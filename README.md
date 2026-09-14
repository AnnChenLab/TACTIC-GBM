# TACTIC-GBM

**Tumor microenvironment–Aware Combination Therapy Index and Characterization for Glioblastoma**

An interactive R Shiny application for exploring candidate genes, single-drug and
drug-combination ("doublet") therapies, and single-cell expression patterns in glioblastoma
(GBM), integrating bulk and single-cell transcriptomic data with survival outcomes.

> TODO: one or two sentences on the scientific motivation / what question this tool helps answer,
> and a link to the associated publication once available.

## Features

- **Single Gene Analysis** — survival stratification by single-gene expression across cohorts
- **Drug-Gene** — gene targets for a selected drug across multiple annotation databases
  (DGIdb, DrugBank, DsigDB, DREIMT), with Euler-diagram overlap visualization
- **Single Drug Analysis** — PCA-based survival stratification for single drugs
- **Top Doublets** — candidate drug-combination ("doublet") ranking table
- **PC1 and Survival** — doublet PC1-based survival stratification
- **Single Cell Heatmap** — single-cell expression heatmaps by cell type/group

## Data sources

Built on de-identified, publicly available cohort data: **TCGA**, **CGGA**, and **ORIEN**.
No patient-identifying information (names, MRNs, dates of birth) is included — only
de-identified research fields (e.g. patient barcode, age, sex, survival/MGMT status).

## Running locally

Requires R (developed against 4.5.3) and the packages listed in the `library()` calls at the
top of `app.R`. From an R console in this directory:

```r
install.packages(c(
  "shiny", "plotly", "ggplot2", "data.table", "dplyr", "coin", "openxlsx", "tidyr", "DT",
  "htmltools", "VennDiagram", "readr", "eulerr", "bslib", "survival", "survminer",
  "gridExtra", "pheatmap", "tibble", "stringr"
))
shiny::runApp()
```

## Running with Docker/Podman

A `Containerfile` is included, targeting `rocker/r-ver:4.5.3` for a reproducible build
environment (system libraries, a Rust toolchain for the `eulerr` package, and a BLAS/LAPACK
pin that avoids a known OpenBLAS illegal-instruction crash on some hosts):

```bash
docker build -t tactic-gbm -f Containerfile .
docker run -p 3838:3838 tactic-gbm
```
Then open `http://localhost:3838`.

## Citation

> TODO: add citation (paper title, authors, journal, DOI) once published.

## License

MIT — see [LICENSE](LICENSE).

## Contact

> TODO: corresponding author / maintainer contact info.
