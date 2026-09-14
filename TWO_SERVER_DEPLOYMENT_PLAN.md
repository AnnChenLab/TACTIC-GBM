# Two-Server Deployment Plan — Rshiny_GBM_dev (TACTIC-GBM)

Companion to `DEPLOYMENT_GUIDE.md` (which documents the current, working, single-server setup on chenlab). This file is the plan to extend that to a second server (chenlab2) and to make both publicly reachable over HTTPS. Nothing in this plan has been executed yet — it's a runbook to work through by hand on each server. Fold the results back into `DEPLOYMENT_GUIDE.md` once done (Phase 5) so there's one durable reference going forward.

**Topology decision**: chenlab2 = staging (test new builds first), chenlab = primary/production (promote verified builds here — it already has a proven, documented process).

**Public-exposure decision**: put a reverse proxy (nginx or Caddy) with TLS in front of the container on both servers, rather than exposing port 3838 directly — this is what actually resolves the "how does this get exposed externally" question that `DEPLOYMENT_GUIDE.md` already flags as unresolved.

---

## Phase 0 — Reconcile inputs before touching either server

- [ ] `git pull` in this folder for the author's (Yuxin's) latest upstream changes.
- [ ] Diff `app.R` after pulling for any new `library(...)` calls not yet in `Containerfile`'s `install.packages()` list — add any missing packages before rebuilding.
- [ ] Resolve the R version mismatch: `Containerfile` pins `rocker/r-ver:4.5.3`; this Mac's installed R is 4.6.0. Confirm with Yuxin whether 4.5.3 is intentionally what's tested against, or bump the base image.
- [ ] Get a one-line confirmation from whoever assembled `data/` that nothing beyond the public TCGA/CGGA/ORIEN-derived fields was merged in. (Columns inspected so far — `patient_barcode`, `age`, `sex`, survival/MGMT fields — look like standard de-identified research data, no MRNs/DOB/names found, but this is about to go publicly live so a human sanity check is cheap insurance.)

## Phase 1 — Per-server facts table

Fill in the right column on first login to chenlab2; leave the chenlab column as-is (already proven, don't touch what works).

| Fact | chenlab (production) | chenlab2 (staging) |
|---|---|---|
| Host | `chenlab.chpc.utah.edu` | **[CONFIRM ON FIRST LOGIN]** |
| Deploy path | `/uufs/chpc.utah.edu/common/home/ychen-group1/shiny_apps/Rshiny_GBM_dev/` | **[CONFIRM]** |
| Service account | shared `podman` account | **[CONFIRM — may not exist; may build under own account instead]** |
| Runtime | podman, rootless | **[CONFIRM podman is installed; fall back to docker if not]** |

## Phase 2 — Staging build & verification on chenlab2 (test first)

1. Fill in the chenlab2 row above before running anything.
2. Sync the project over (same exclude list as the existing chenlab rsync command in `DEPLOYMENT_GUIDE.md`):
   ```bash
   rsync -avzP --exclude '.git' --exclude '.Rproj.user' --exclude '.DS_Store' \
     --exclude '.RData' --exclude '.Rhistory' --exclude 'shiny_run.log' \
     "/Users/chloetran/Library/CloudStorage/Box-Box/shiny_app/Rshiny_GBM_dev/" \
     <user>@<chenlab2-host>:<chenlab2-deploy-path>/
   ```
3. Build a date-tagged image (keep the date-tag-not-`latest` convention so rollback stays one command):
   ```bash
   podman build -t gbm-shiny:$(date +%Y%m%d) -f Containerfile .
   ```
4. Run it **bound to localhost only** first — not yet public:
   ```bash
   podman run -d --name gbm-shiny-staging --restart=always -p 127.0.0.1:3838:3838 gbm-shiny:<tag>
   ```
5. Verify before trusting the build:
   ```bash
   curl localhost:3838                                   # from the node itself
   ssh -L 3838:localhost:3838 <user>@<chenlab2-host>      # then open http://localhost:3838 locally
   ```
   Click through the app — data loads, plots render, tables populate.
6. Only after that passes, move to Phase 3. (chenlab2 can stay internal-only during staging, or also go public via the same proxy config if it needs to be reachable by others during testing — your call at execution time.)

## Phase 3 — Reverse proxy + TLS (same pattern on both servers)

- Prefer **Caddy** for automatic Let's Encrypt TLS if the server has a public DNS name and ports 80/443 open. Fall back to **nginx + certbot** if Caddy isn't already standard/approved on CHPC.
- The proxy config **must** forward WebSocket upgrade headers to `127.0.0.1:3838`:
  ```
  proxy_set_header Upgrade $http_upgrade;
  proxy_set_header Connection "upgrade";
  ```
  (Caddy does this automatically; for nginx it must be set explicitly.) Shiny's reactivity runs over a websocket — without this, the page loads but clicking anything does nothing, which is the single most common Shiny-behind-proxy failure.
- In parallel (not blocking the technical work): file a request with CHPC support to confirm inbound 80/443 is actually permitted on chenlab/chenlab2's network policy — `DEPLOYMENT_GUIDE.md` already flags this as unconfirmed even for chenlab today.
- Confirm whether CHPC expects the proxy to run as another rootless container (a podman pod) or as a node-level system service, and match whatever convention is already standard there rather than improvising.

## Phase 4 — Promote the verified build to chenlab (production)

1. Repeat the same rsync + `podman build` with the same date tag **on chenlab** (rebuild from source rather than copying the image artifact between servers — keeps the existing tag/rollback model intact).
2. Follow `DEPLOYMENT_GUIDE.md`'s existing "Updating the server to a new version" steps verbatim — that part already works, no changes needed:
   ```bash
   podman stop gbm-shiny 2>/dev/null; podman rm gbm-shiny 2>/dev/null
   podman run -d --name gbm-shiny --restart=always -p 3838:3838 gbm-shiny:<new-tag>
   ```
3. Apply the same reverse-proxy/TLS config from Phase 3.
4. Smoke-test chenlab's public URL from **off-node** (actual public internet, not just from within CHPC) to confirm real public reachability.

## Phase 5 — Fold this back into DEPLOYMENT_GUIDE.md

Once the above is done and working, merge into `DEPLOYMENT_GUIDE.md`: the per-server facts table (Phase 1), the staging→production promotion workflow (Phases 2 & 4), and the reverse-proxy setup + WebSocket gotcha (Phase 3). Then this plan file can be deleted — `DEPLOYMENT_GUIDE.md` becomes the one file a future maintainer needs.

---

## Open items this plan can't resolve without server access

- chenlab2's real hostname, deploy path, and container runtime — must be filled in on first login.
- Whether CHPC's network policy permits inbound 80/443 — a support ticket, not resolvable in advance.
- R version reconciliation (4.5.3 pinned vs. 4.6.0 local dev) — needs Yuxin's input.
- No authentication layer exists in the app. Since the underlying data is de-identified public research data (TCGA/CGGA/ORIEN), this plan treats that as acceptable for "public accessing" — flagged here in case that assumption is wrong.

## Verification checklist

- [ ] Staging (chenlab2): container starts; `curl localhost:3838` returns 200; SSH-tunneled browser test loads the app and its plots/tables render.
- [ ] Reverse proxy: app stays *interactive* after load (a control click actually reacts), not just that the HTML page loads.
- [ ] Production (chenlab): same smoke test, run from outside the CHPC network to confirm actual public reachability.
- [ ] Rollback safety: previous date-tagged image still present after each promotion (`podman images`), so reverting is a one-line `podman run` to the prior tag.
