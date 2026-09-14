# Deployment Guide — Rshiny_GBM_dev (TACTIC-GBM)

This file documents the current, working, single-server setup on `chenlab` (production). A
companion file, `TWO_SERVER_DEPLOYMENT_PLAN.md`, is the in-progress (not yet executed) plan to
extend this to a second staging server (`chenlab2`) plus a reverse proxy/TLS layer for public
access — read that file for the next-steps roadmap; this file is the ground truth for how the
existing deployment actually works today.

## What this project is

An R Shiny app (glioblastoma drug-combination explorer) maintained by the author at
`https://github.com/yuxin6058887/TACTIC-GBM.git`. This local folder (`Rshiny_GBM_dev/`) is a
git checkout of that repo, with two extra local-only files layered on top that are **not**
part of the upstream repo:

| File | Tracked by upstream git? | Purpose |
|---|---|---|
| `app.R` | Yes | App entrypoint (ui + server + `shinyApp()`) |
| `data/` | Yes | ~88MB of `.rds`/`.csv`/`.xlsx` data files, read via relative paths |
| `Containerfile` | No | Our own Podman/Docker build definition |
| `.containerignore` | No | Excludes local cruft from the build context |
| `DEPLOYMENT_GUIDE.md` | No | This file |

Future updates from the author: `cd Rshiny_GBM_dev && git pull`, then rebuild the image (see below).
Only `Containerfile`/`.containerignore` need manual edits, and only if the author adds new
R package dependencies (new `library()` calls) — check the diff after pulling.

## Containerfile summary

- Base image: `rocker/r-ver:4.5.3` (matches the R version the app was developed against; upstream
  doesn't pin one, so this is our own choice — revisit if the author specifies otherwise)
- System packages: `curl cmake libuv1-dev build-essential gfortran libcurl4-openssl-dev
  libssl-dev libxml2-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev
  libjpeg-dev libicu-dev libcairo2-dev zlib1g-dev`
- Rust via `rustup` (the `eulerr` package needs rustc ≥1.88; Ubuntu's apt `rustc` is older)
- R packages installed directly via `install.packages()` (no `renv.lock` — dropped since it
  wasn't used by the build): `shiny plotly ggplot2 data.table dplyr coin openxlsx tidyr DT
  htmltools VennDiagram readr eulerr bslib survival survminer gridExtra pheatmap tibble
  stringr` (+ ~110 transitive dependencies resolved fresh at build time)
- Entrypoint: `shiny::runApp('/app/app.R', host='0.0.0.0', port=3838)`, port `3838` exposed

## Gotchas hit during setup (for future reference)

1. **`eulerr` needs a modern Rust toolchain.** Fixed by installing via `rustup` instead of
   `apt install rustc`.
2. **Local Mac testing must build natively (arm64), not `--platform linux/amd64`.** Forcing the
   amd64 platform on Apple Silicon runs the whole build under QEMU emulation, and `rustc`
   segfaults under that emulation (`qemu: uncaught target signal 11`). This only affects local
   testing — the real CHPC build runs on native x86_64 hardware and won't hit this.
   `podman build --platform linux/arm64 -t gbm-shiny:local -f Containerfile .`
3. **Box Drive's virtiofs sharing can corrupt file reads mid-build** (`input/output error`,
   `unexpected EOF` on `COPY app.R`). If a local test build fails this way, copy the project to
   plain local disk first and build from there instead of directly off the Box-synced folder:
   ```bash
   rsync -a --exclude '.git' "$PWD/" /tmp/Rshiny_GBM_dev_build/
   cd /tmp/Rshiny_GBM_dev_build && podman build --platform linux/arm64 -t gbm-shiny:local -f Containerfile .
   ```
   This is a local-Mac-only issue — CHPC deployment goes over SSH to a real filesystem, not Box Drive.
4. **OpenBLAS illegal-instruction crash (exit 132 / SIGILL) inside `grDevices::convertColor`,
   triggered by `eulerr`'s Euler-diagram rendering (Venn/Euler tab).** OpenBLAS's runtime
   CPU-feature detection can pick an instruction set the actual host doesn't support, crashing
   any BLAS/LAPACK-backed call — this hit in production on a real server, not just locally. Fixed
   by pinning BLAS/LAPACK to the reference implementation, now baked into `Containerfile` right
   after the package install step:
   ```dockerfile
   RUN update-alternatives --set libblas.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/blas/libblas.so.3 && \
       update-alternatives --set liblapack.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3
   ```
   **This bit us once already by being applied only to a server's copy of `Containerfile` and not
   carried back into this local checkout** — the next rsync from this Mac silently overwrote the
   fix. It's now in the source-of-truth `Containerfile` here, so it survives future syncs. If a
   future rebuild ever crashes with exit 132 again, this is the first thing to check.

## Server-side facts (CHPC)

- Host: `chenlab.chpc.utah.edu`, user `u6057318`
- Deployment path: `/uufs/chpc.utah.edu/common/home/ychen-group1/shiny_apps/Rshiny_GBM_dev/`
- Builds run under the shared `podman` service account (`sudo su - podman`) — this account holds
  the subuid/subgid mapping CHPC's shared filesystem needs for rootless container builds
- Container listens on port `3838` internally, but is published on **host port `3839`**
  (`-p 3839:3838`) — confirmed from the actual running container on chenlab (`podman ps -a`
  showed `0.0.0.0:3839->3838/tcp`). Use `3839` for any host-side command (`curl localhost:3839`,
  firewall rules, reverse-proxy upstream target); `3838` only refers to the container's internal
  port. How `3839` itself gets exposed off-node (reverse proxy, firewall rule) is still
  unconfirmed with CHPC support/docs.
- Rootless-Podman gotcha: files touched by the `podman` account can end up with different
  UID/GID ownership. Treat the account's copy as build-only staging, keep source-of-truth edits
  under your own `u6057318` account.

## Updating the server to a new version

```bash
# 1. Sync the updated project to CHPC (run from this Mac)
rsync -avzP --exclude '.git' --exclude '.Rproj.user' --exclude '.DS_Store' \
  --exclude '.RData' --exclude '.Rhistory' --exclude 'shiny_run.log' \
  "/Users/chloetran/Library/CloudStorage/Box-Box/shiny_app/Rshiny_GBM_dev/" \
  u6057318@chenlab.chpc.utah.edu:/uufs/chpc.utah.edu/common/home/ychen-group1/shiny_apps/Rshiny_GBM_dev/

# 2. On the server, as the podman service account, rebuild with a new date-tagged image
sudo su - podman
cd /uufs/chpc.utah.edu/common/home/ychen-group1/shiny_apps/Rshiny_GBM_dev
podman build -t gbm-shiny:$(date +%Y%m%d) -f Containerfile .
# (tag by date, not `latest`, so you can roll back to the previous tag if something breaks)

# 3. Swap the running container over to the new image
podman stop gbm-shiny 2>/dev/null; podman rm gbm-shiny 2>/dev/null
podman run -d --name gbm-shiny --restart=always -p 3839:3838 gbm-shiny:<new-tag>
```

## Operating the running container

```bash
podman logs -f gbm-shiny            # tail logs
podman logs --tail 200 gbm-shiny    # recent logs only
podman restart gbm-shiny            # restart without rebuilding
curl localhost:3839                 # sanity check from within the node (host port, not 3838)
```

## Local testing (before pushing to CHPC)

```bash
cd Rshiny_GBM_dev
podman build --platform linux/arm64 -t gbm-shiny:local -f Containerfile .
podman run --rm -d --name gbm-shiny-test -p 3838:3838 gbm-shiny:local
# wait ~15-20s for startup (loads several large .rds files), then open http://localhost:3838
podman stop gbm-shiny-test
```
