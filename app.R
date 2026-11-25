library(shiny)
library(plotly)
library(ggplot2)
library(data.table)
library(dplyr)
library(coin)
library(openxlsx)
library(tidyr)
library(DT)
library(htmltools)
library(VennDiagram)
library(grid)
library(readr)
library(eulerr)

####### Here is the part for Shiny UI #######
library(shiny)
library(bslib)
library(survival)
library(survminer)
library(dplyr)
library(gridExtra)
library(pheatmap)
library(tibble)
library(stringr)

####### Contact form settings #######
web3forms_access_key <- "0f471c59-725c-41fa-bb79-660c7b6c7461"   



tcga_clinical_subtype = readRDS("data/tcga_clinical_subtype.rds")
cgga325_clinical_subtype= readRDS("data/CGGA_clinical_325_subtype.rds")
#cgga693_clinical_subtype = readRDS("data/CGGA_clinical_693_subtype.rds")


gbm_tcga_seq_final_geneanlysis <- readRDS("data/gbm_tcga_seq_final_geneanlysis_new.rds")
#CGGA_rna_seq_final_693_geneanlysis <- readRDS("data/CGGA_693_readcount_Dseq2.rds")
CGGA_rna_seq_final_325_geneanlysis <- readRDS("data/CGGA_325_readcount_Dseq2.rds")

####### Primary dataset for use ###########
gbm_tcga_seq_final_geneanlysis_primary = gbm_tcga_seq_final_geneanlysis %>% filter(patient_barcode %in% tcga_clinical_subtype$patient_barcode) 
#CGGA_rna_seq_final_693_geneanlysis_primary = CGGA_rna_seq_final_693_geneanlysis %>% filter(patient_barcode %in% cgga693_clinical_subtype$patient_barcode)
CGGA_rna_seq_final_325_geneanlysis_primary = CGGA_rna_seq_final_325_geneanlysis %>% filter(patient_barcode %in% cgga325_clinical_subtype$patient_barcode)

gbm_datasets <- c(
  "GSE103224",
  "GSE131928_10X",
  "GSE131928_Smartseq2",
  "GSE138794",
  "GSE139448",
  "GSE141982",
  "GSE148842",
  "GSE141383",
  "GSE141460",
  "GSE162631",
  "GSE163108_10X",
  "GSE163108_Smartseq2"
)

###### Final result listed here ##########
Doublet_all_with_coef = read.xlsx("data/Doublet_all_with_coef.xlsx")
candidate_doublet = read_csv("data/s_table2.csv")

obj_names <- load("data/drug_db_ls_v4.RData")
drug_ls = c()
for (i in names(drug_db_ls_new)) {
  drug = colnames(drug_db_ls_new[[i]])
  drug_ls = c(drug_ls,drug)
  drug_ls = drug_ls[drug_ls != "gene"]
}

drug_set_ls <- lapply(drug_ls, function(drug) {
  names(drug_db_ls_new)[
    sapply(drug_db_ls_new, function(x) drug %in% colnames(x))
  ]
})
names(drug_set_ls) <- drug_ls

drug_set_ls_grouped <- lapply(
  drug_set_ls,
  function(x) {
    x_grouped <- ifelse(
      grepl("dsigdb", x, ignore.case = TRUE),
      "DsigDB",
      x
    )
    
    unique(x_grouped)
  }
)
names(drug_set_ls_grouped) <- drug_ls


unique_genes <- Doublet_all_with_coef%>%
  select(target_gene_drug1, target_gene_drug2) %>%
  pivot_longer(cols = everything()) %>%
  separate_rows(value, sep = ",\\s*") %>%
  distinct(value) %>%
  pull(value)


tcga_gene = readRDS("data/gbm_tcga_seq_final_geneanlysis.rds")
gene_ls = colnames(tcga_gene)
drug_gene_ls = readRDS("data/drug_gene_list.rds")
drug_list = readRDS("data/drug_list.rds")

drug_analysis_unique_df_cluster_combined = read.csv("data/single_drug_result_unique_cluster_1020.csv")
drug_seed_df = drug_analysis_unique_df_cluster_combined %>%
  filter(p_tcga<0.05 & p_cgga325<0.05)
exist_drugs <- c("Temozolomide", "Bevacizumab", "Lomustine", "Irinotecan", "Pembrolizumab")
matched_drugs <- unique(drug_list[str_detect(drug_list,regex(str_c(exist_drugs, collapse = "|"), ignore_case = TRUE))])

drug_seed_ls = c(drug_seed_df$drug,matched_drugs,"radio_rsi")

sorted_gene_sets <- lapply(drug_gene_ls, function(genes) sort(genes))
gene_set_keys <- sapply(sorted_gene_sets, paste, collapse = ",")
unique_indices <- !duplicated(gene_set_keys)
unique_drug_genes <- drug_gene_ls[unique_indices]
unique_drug_genes[["radio_rsi"]] <- c("AR", "JUN", "STAT1", "PRKCB", "RELA", "ABL1", "SUMO1", "PAK2", "HDAC1", "IRF1")

drug_gene_tab_choices <- sort(
  unique(
    unlist(
      lapply(drug_db_ls_new, function(dat) {
        setdiff(colnames(dat), "gene")
      }),
      use.names = FALSE
    )
  )
)

nice_surv_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

merge_by_group = read.csv("data/heatmap_by_group0323.csv")
merge_by_type = read.csv("data/heatmap_by_type0323.csv")
merge_by_group <- merge_by_group %>%
  mutate(
    dataset = str_extract(
      file_name,
      "GSE\\d+(?:_(?:10X|Smartseq2))?"
    )
  )
merge_by_type <- merge_by_type %>%
  mutate(
    dataset = str_extract(
      file_name,
      "GSE\\d+(?:_(?:10X|Smartseq2))?"
    )
  )
merge_by_group_5 = merge_by_group %>% select(-Fibroblast)

merge_by_group_gene = read_rds("data/heatmap_gene_by_group.rds")
merge_by_type_gene = read_rds("data/heatmap_gene_by_type.rds")

merge_by_group_gene <- merge_by_group_gene %>%
  mutate(
    dataset = str_extract(file_name, "GSE\\d+(?:_(?:10X|Smartseq2))?")
  )

merge_by_type_gene <- merge_by_type_gene %>%
  mutate(
    dataset = str_extract(file_name, "GSE\\d+(?:_(?:10X|Smartseq2))?")
  )

####### Manuscript-featured doublets ###########
# The three doublets discussed in detail in the manuscript Results
# (Table 1 ranks 1, 5 and 4). Names are matched to the app's data
# case- and punctuation-insensitively, so "Gw-559090" matches "GW-559090".
featured_doublets <- tibble::tribble(
  ~label,                        ~drug1,       ~drug2,                 ~table1_rank, ~rra_p,
  "GW-559090 + CHEMBL1321401",   "Gw-559090",  "Chembl:Chembl1321401", 1,            "<0.001",
  "PD-184352 + SAG1.3",          "Pd-184352",  "Sag1.3",               5,            "0.088",
  "RWJ-67657 + AS-1409",         "Rwj-67657",  "As-1409",              4,            "0.088"
)

normalize_drug_name <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))

# Return the spelling used in `choices` for each name in `wanted` (NA if absent)
match_drug_names <- function(wanted, choices) {
  choices <- as.character(choices)
  vapply(wanted, function(w) {
    hit <- choices[normalize_drug_name(choices) == normalize_drug_name(w)]
    if (length(hit)) hit[1] else NA_character_
  }, character(1), USE.NAMES = FALSE)
}

# Put featured drugs first, as their own option group in the dropdown
prioritize_choices <- function(choices, featured,
                               featured_label = "Discussed in manuscript",
                               rest_label = "All drugs") {
  choices <- unique(as.character(choices))
  choices <- choices[!is.na(choices)]
  hits <- unique(stats::na.omit(match_drug_names(featured, choices)))
  if (length(hits) == 0) return(choices)
  out <- list(as.list(hits), as.list(setdiff(choices, hits)))
  names(out) <- c(featured_label, rest_label)
  out
}

featured_all_drugs <- unique(c(featured_doublets$drug1, featured_doublets$drug2))

####### Tutorial screenshots ###########
# Save four screenshots in ONE folder and put its FULL path below, e.g.
#   "data/tutorial_screenshots"
# Expected names (.png, .jpg or .jpeg; upper/lower case doesn't matter):
#   step1_top_combination, step1_evidence, step2_doublet, step2_singledrug,
#   step2_singlegene, step3_sc_bydrug, step3_sc_bygene
# The pictures are embedded directly in the page, so no www/ folder is needed.
# Restart the app after adding or replacing a picture.
tutorial_shot_dir <- "data/tutorial_screenshots"   # relative to the app folder

# Find the screenshot folder: the path above as given, then next to the app
# file, then in the working directory.
resolve_tutorial_dir <- function(dir) {
  app_dir <- tryCatch({
    f <- NULL
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)$ofile
      if (!is.null(f)) break
    }
    if (is.null(f) && requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      f <- rstudioapi::getSourceEditorContext()$path
    }
    if (!is.null(f) && nzchar(f)) dirname(normalizePath(f)) else NULL
  }, error = function(e) NULL)
  
  candidates <- unique(c(
    path.expand(dir),
    if (!is.null(app_dir)) file.path(app_dir, basename(dir)),
    file.path(getwd(), basename(dir))
  ))
  hit <- candidates[dir.exists(candidates)]
  if (length(hit)) normalizePath(hit[1]) else normalizePath(path.expand(dir), mustWork = FALSE)
}
tutorial_shot_dir <- resolve_tutorial_dir(tutorial_shot_dir)

# Match by name, ignoring case, spaces/punctuation and the extension, so
# "Step3_SingleCell.PNG" or "step3_singlecell (1).png" still match.
find_tutorial_shot <- function(file) {
  if (!dir.exists(tutorial_shot_dir)) return(NA_character_)
  key <- function(x) gsub("[^a-z0-9]", "", tolower(tools::file_path_sans_ext(basename(x))))
  files <- list.files(tutorial_shot_dir, full.names = TRUE)
  files <- files[tolower(tools::file_ext(files)) %in% c("png", "jpg", "jpeg")]
  hit <- files[startsWith(key(files), key(file))]
  if (length(hit)) hit[1] else NA_character_
}

local({
  shots <- c("step1_top_combination", "step1_evidence", "step2_doublet",
             "step2_singledrug", "step2_singlegene", "step3_sc_bydrug", "step3_sc_bygene")
  found <- !is.na(vapply(shots, find_tutorial_shot, character(1)))
  message("Tutorial screenshots: looking in ", tutorial_shot_dir,
          if (!dir.exists(tutorial_shot_dir)) "  (folder NOT found)" else "",
          "\n  found: ", if (any(found)) paste(shots[found], collapse = ", ") else "none",
          if (any(!found)) paste0("\n  missing: ", paste(shots[!found], collapse = ", ")) else "")
})

tutorial_shot <- function(file, caption) {
  path <- find_tutorial_shot(file)
  if (!is.na(path)) {
    mime <- if (tolower(tools::file_ext(path)) == "png") "image/png" else "image/jpeg"
    tags$figure(
      class = "tutorial-shot",
      tags$img(src = base64enc::dataURI(file = path, mime = mime), alt = caption,
               title = "Click to enlarge", onclick = "tutorialZoom(this)"),
      tags$figcaption(caption)
    )
  } else {
    div(class = "tutorial-shot tutorial-shot-missing",
        icon("image"), " Screenshot not found: ", tags$code(file),
        div(style = "font-size: 12px; margin-top: 4px;",
            "Looking in ", tags$code(tutorial_shot_dir),
            if (!dir.exists(tutorial_shot_dir)) {
              " (this folder does not exist; set tutorial_shot_dir to the full path of your screenshot folder)"
            } else {
              shown <- list.files(tutorial_shot_dir)
              tagList(br(), "Files in that folder: ",
                      if (length(shown)) paste(shown, collapse = ", ") else "(empty)")
            }))
  }
}

# Several screenshots side by side (wraps on narrow screens)
tutorial_gallery <- function(files, captions) {
  div(class = "tutorial-gallery", Map(tutorial_shot, files, captions))
}

# Selection boxes: clicking clears the current name (shown grey as a hint)
# so users can type straight away; clicking away without choosing restores it.
# The clear is silent, so the server keeps the old value until a new pick.
click_to_type_options <- function(placeholder = "Type to search...") {
  list(
    placeholder = placeholder,
    onInitialize = I("function() {
      var self = this;
      self.on('focus', function() {
        var v = self.getValue();
        if (!v) return;
        var opt = self.options[v];
        self._prevValue = v;
        self.settings.placeholder = opt ? opt[self.settings.labelField] : v;
        self.clear(true);
        self.updatePlaceholder();
      });
      self.on('blur', function() {
        if (!self.getValue() && self._prevValue) {
          self.setValue(self._prevValue, true);
        }
        self._prevValue = null;
      });
    }")
  )
}

doublet_drug1_choices <- prioritize_choices(drug_seed_ls, featured_doublets$drug1)
doublet_drug2_choices <- prioritize_choices(unique(Doublet_all_with_coef$Drug2), featured_doublets$drug2)
single_drug_choices   <- prioritize_choices(drug_list, featured_all_drugs)
heatmap_drug_choices  <- prioritize_choices(unique(merge_by_group$drug), featured_all_drugs)

# Resolve featured names to the spellings actually present in each selector
featured_doublets <- featured_doublets %>%
  mutate(
    app_drug1 = match_drug_names(drug1, drug_seed_ls),
    app_drug2 = match_drug_names(drug2, unique(Doublet_all_with_coef$Drug2))
  )
local({
  missing <- c(featured_doublets$drug1[is.na(featured_doublets$app_drug1)],
               featured_doublets$drug2[is.na(featured_doublets$app_drug2)])
  if (length(missing)) {
    warning("Featured drug(s) not found in the Doublet Survival selectors: ",
            paste(missing, collapse = ", "))
  }
})

####### Score definitions (single source for tooltips + tutorial) ###########
score_definitions <- list(
  `Drug 1` = "Seed therapy: one of the 48 single agents (29 with survival association in both cohorts, plus standard-of-care GBM agents and radiotherapy RSI) paired with every other drug.",
  `Drug 2` = "Partner drug paired with the seed therapy.",
  `Drug 1 Source` = "Database the Drug 1 target genes come from (DGIdb, DrugBank, DREIMT, or the DSigDB signature set dsigdb_d3).",
  `Drug 2 Source` = "Database the Drug 2 target genes come from (DGIdb, DrugBank, DREIMT, or the DSigDB signature set dsigdb_d3).",
  `Fisher-combined P-value` = "Multivariable Cox model with PC1 of Drug 1 and Drug 2 target genes, adjusted for age and GBM subtype, fit separately in TCGA and CGGA325. The two cohort p-values are combined with Fisher's product method. Lower = stronger. This is the 'Fisher Product P-value' column in manuscript Table 1.",
  `Clinical Association Score` = "-log10(Fisher-combined P-value), as defined in the manuscript. Higher = stronger combined survival association. 1.30 corresponds to P = 0.05, 2 to P = 0.01, 3 to P = 0.001.",
  `TME Specificity Score` = "Sum of the two drugs' specificity scores. Each drug's score is its maximum target score across single-cell populations; higher = more cell-type-specific targeting.",
  `Literature Semantic Score` = "Sum of the two drugs' literature scores. Each drug's score is the mean PubMedBERT similarity plus a disease-specificity bonus (GBM > glioma/brain tumor > cancer) over its top 5 PubMed records. Higher = more literature support.",
  `Reciprocal Rank Fusion Score` = "Reciprocal Rank Fusion: sum of 1/(k + rank) over the three evidence rankings. Higher = better. Shown as a sensitivity check; the manuscript ranks doublets by RRA.",
  `RRA P-value (TACTIC Score)` = "Robust Rank Aggregation across the clinical, TME and literature rankings. Probability-like: lower = more consistently near the top of all three lists. The table is sorted by this column, as in manuscript Table 1."
)

####### Clinical Association Score: align with manuscript Table 1 ###########
# Manuscript: Clinical association score = -log10(Fisher-combined P).
# The stored column Fisher_combined_tcga_cgga325 may be on another scale
# (raw P, -ln P, or the Fisher chi-square statistic), which shows up as a
# constant-factor difference from Table 1. We calibrate against Table 1.
# Set to "p", "neg_log10_p", "neg_ln_p" or "fisher_chisq_df4" to override.
clinical_score_scale <- "auto"

manuscript_table1 <- tibble::tribble(
  ~drug1,                     ~drug2,                  ~clinical_score,
  "Gw-559090",                "Chembl:Chembl1321401",  3.81,
  "2-Aminoethanimidic Acid",  "Chembl:Chembl1321401",  1.86,
  "Sag1.3",                   "E-6201",                2.65,
  "Rwj-67657",                "As-1409",               2.45,
  "Pd-184352",                "Sag1.3",                2.97,
  "Sag1.3",                   "Apatorsen",             2.90,
  "Gw-559090",                "Chembl:Chembl579882",   1.30,
  "Rwj-67657",                "Dihydrotanshinone I",   4.17,
  "Rwj-67657",                "Methoxsalen_Pc3_Down",  2.62,
  "Rwj-67657",                "Nabumetone_Pc3_Down",   3.77
)

clinical_score_converters <- list(
  neg_log10_p      = function(v) v,
  p                = function(v) -log10(v),
  neg_ln_p         = function(v) v / log(10),
  fisher_chisq_df4 = function(v) -pchisq(v, df = 4, lower.tail = FALSE, log.p = TRUE) / log(10)
)

detect_clinical_scale <- function(values, drug1, drug2) {
  pair_key <- function(a, b) {
    a <- normalize_drug_name(a); b <- normalize_drug_name(b)
    paste(pmin(a, b), pmax(a, b))
  }
  idx <- match(pair_key(manuscript_table1$drug1, manuscript_table1$drug2),
               pair_key(drug1, drug2))
  ok <- !is.na(idx)
  
  if (!any(ok)) {
    guess <- if (all(values > 0 & values <= 1, na.rm = TRUE)) "p" else "neg_log10_p"
    warning("No Table 1 doublets found in candidate table; assuming Clinical column is '", guess, "'.")
    return(guess)
  }
  
  v <- values[idx[ok]]
  target <- manuscript_table1$clinical_score[ok]
  err <- vapply(names(clinical_score_converters), function(nm) {
    if (nm == "p" && any(v <= 0 | v > 1, na.rm = TRUE)) return(Inf)
    s <- suppressWarnings(clinical_score_converters[[nm]](v))
    stats::median(abs(s - target), na.rm = TRUE)
  }, numeric(1))
  
  best <- names(which.min(err))
  message(sprintf(
    "Clinical Association Score: stored column looks like '%s' (median |diff| vs Table 1 = %.3f on %d doublets).",
    best, err[[best]], sum(ok)
  ))
  if (err[[best]] > 0.05) {
    warning("Clinical Association Score still differs from manuscript Table 1 after conversion; check the source column.")
  }
  best
}

clinical_scale_used <- if (identical(clinical_score_scale, "auto")) {
  detect_clinical_scale(candidate_doublet$Fisher_combined_tcga_cgga325,
                        candidate_doublet$Drug1, candidate_doublet$Drug2)
} else {
  clinical_score_scale
}

candidate_doublet <- candidate_doublet %>%
  mutate(
    clinical_score = clinical_score_converters[[clinical_scale_used]](Fisher_combined_tcga_cgga325),
    clinical_p     = 10^(-clinical_score)
  )

# Doublets in the same order as the Top Combinations table (by RRA)
candidate_ranked <- candidate_doublet %>%
  arrange(RRA_score) %>%
  mutate(
    rra_rank   = row_number(),
    pair_label = paste0(rra_rank, ". ", Drug1, " + ", Drug2)
  )
doublet_pair_choices <- candidate_ranked$pair_label

####### PubMed IDs by drug ###########
pubmed_path <- "data/Pubmed Retrieve ID by Drug.xlsx"

pubmed_by_drug <- if (file.exists(pubmed_path)) {
  openxlsx::read.xlsx(pubmed_path)
} else {
  warning("PubMed file not found: ", pubmed_path)
  data.frame(drug = character(0), targets = character(0), pmid_list = character(0),
             drug_name_pmids = character(0), target_gene_pmids = character(0),
             top5_pmid_list = character(0))
}

split_ids <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
  v <- trimws(strsplit(as.character(x), ";")[[1]])
  unique(v[nzchar(v)])
}

# One row per drug x PMID (dominant cell group is added further down)
pmid_long <- dplyr::bind_rows(lapply(seq_len(nrow(pubmed_by_drug)), function(i) {
  r   <- pubmed_by_drug[i, ]
  ids <- split_ids(r$pmid_list)
  if (!length(ids)) return(NULL)
  tibble::tibble(
    Drug           = r$drug,
    `Target genes` = r$targets,
    PMID           = ids,
    `PubMed link`  = paste0("https://pubmed.ncbi.nlm.nih.gov/", ids, "/")
  )
}))
if (nrow(pmid_long) == 0) {
  pmid_long <- tibble::tibble(Drug = character(0), `Target genes` = character(0),
                              PMID = character(0), `PubMed link` = character(0))
}

pmid_drug_choices <- prioritize_choices(sort(unique(pubmed_by_drug$drug)), featured_all_drugs)

# Row of a data frame whose `col` matches `drug` (case/punctuation-insensitive)
row_for_drug <- function(df, col, drug) {
  if (!nrow(df)) return(df)
  df[normalize_drug_name(df[[col]]) == normalize_drug_name(drug), , drop = FALSE][1, , drop = FALSE]
}

####### Dominant cell group per drug (for the evidence package) ###########
# Mean % of cells expressing the drug's target genes across GBM single-cell
# datasets, per merged cell group; the dominant group is the highest one.
cell_group_cols <- names(merge_by_group_5)[vapply(merge_by_group_5, is.numeric, logical(1))]
cell_group_cols <- setdiff(cell_group_cols, c("X", "...1"))

dominant_by_drug <- if (length(cell_group_cols) == 0) {
  warning("No numeric cell-group columns in merge_by_group_5; dominant cell groups will be blank.")
  tibble::tibble(drug = character(0), dominant_group = character(0))
} else merge_by_group_5 %>%
  filter(dataset %in% gbm_datasets) %>%
  group_by(drug) %>%
  summarise(across(all_of(cell_group_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-drug, names_to = "cell_group", values_to = "pct") %>%
  filter(!is.na(pct)) %>%
  group_by(drug) %>%
  arrange(desc(pct), .by_group = TRUE) %>%
  summarise(
    dominant_group = first(cell_group),
    .groups = "drop"
  )

dominant_group_for <- function(drugs) {
  dominant_by_drug$dominant_group[match(normalize_drug_name(drugs),
                                        normalize_drug_name(dominant_by_drug$drug))]
}

pmid_long <- pmid_long %>%
  mutate(`Dominant cell group` = dominant_group_for(Drug)) %>%
  relocate(`Dominant cell group`, .after = `Target genes`)

target_genes_for <- function(drug) {
  r <- row_for_drug(pubmed_by_drug, "drug", drug)
  if (nrow(r) && !is.na(r$targets[1]) && nzchar(r$targets[1])) return(r$targets[1])
  hit <- names(drug_gene_ls)[normalize_drug_name(names(drug_gene_ls)) == normalize_drug_name(drug)]
  if (length(hit)) paste(sort(unique(unlist(drug_gene_ls[[hit[1]]]))), collapse = "; ") else NA_character_
}

# Evidence package: one row per drug of the selected doublet
build_evidence_package <- function(pair_row) {
  dplyr::bind_rows(lapply(c(1, 2), function(k) {
    drug   <- pair_row[[paste0("Drug", k)]]
    source <- pair_row[[paste0("drug_set_drug", k)]]
    dom    <- row_for_drug(dominant_by_drug, "drug", drug)
    pm     <- row_for_drug(pubmed_by_drug, "drug", drug)
    pick   <- function(df, col) if (nrow(df) && col %in% names(df)) df[[col]][1] else NA
    all_ids <- split_ids(pick(pm, "pmid_list"))
    tibble::tibble(
      Doublet                         = paste(pair_row$Drug1, "+", pair_row$Drug2),
      `RRA rank`                      = pair_row$rra_rank,
      Role                            = paste("Drug", k),
      Drug                            = drug,
      `Target gene source`            = source,
      `Target genes`                  = target_genes_for(drug),
      `Dominant cell group`           = pick(dom, "dominant_group"),
      `Fisher-combined P-value`       = signif(pair_row$clinical_p, 3),
      `Clinical Association Score`    = round(pair_row$clinical_score, 2),
      `TME Specificity Score`         = round(pair_row$TME_score, 3),
      `Literature Semantic Score`     = round(pair_row$literature_score, 3),
      `Reciprocal Rank Fusion Score`  = round(pair_row$RRF_score, 3),
      `RRA P-value (TACTIC Score)`    = signif(pair_row$RRA_score, 3),
      `Key PMIDs (top 5)`             = paste(split_ids(pick(pm, "top5_pmid_list")), collapse = "; "),
      `All PMIDs`                     = paste(all_ids, collapse = "; "),
      `Number of PMIDs`               = length(all_ids)
    )
  }))
}

safe_file_part <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

tutorial_heading <- function(icon_name, text) {
  h5(class = "tutorial-heading", icon(icon_name), text)
}

# Convert original database names into the three displayed source groups
group_database_name <- function(database_name) {
  dplyr::case_when(
    grepl("^DGIdb$", database_name, ignore.case = TRUE) ~ "DGIdb",
    grepl("^DrugBank$", database_name, ignore.case = TRUE) ~ "DrugBank",
    grepl("dsigdb", database_name, ignore.case = TRUE) ~ "DsigDB",
    grepl("dreimt", database_name, ignore.case = TRUE) ~ "DREIMT",
    TRUE ~ NA_character_
  )
}


# Obtain target genes for one drug from one grouped source
get_drug_genes_by_source <- function(drug_name, source_name, database_list) {
  
  relevant_databases <- names(database_list)[
    vapply(
      names(database_list),
      function(database_name) {
        identical(group_database_name(database_name), source_name)
      },
      logical(1)
    )
  ]
  
  genes <- unlist(
    lapply(relevant_databases, function(database_name) {
      
      dat <- database_list[[database_name]]
      
      # Skip databases that do not contain this drug
      if (!drug_name %in% colnames(dat)) {
        return(character(0))
      }
      
      gene_column <- if ("gene" %in% colnames(dat)) {
        as.character(dat$gene)
      } else {
        rownames(dat)
      }
      
      drug_values <- dat[[drug_name]]
      
      # Keep genes linked to the selected drug
      gene_column[
        !is.na(drug_values) &
          drug_values == 1
      ]
    }),
    use.names = FALSE
  )
  
  unique(na.omit(genes))
}

wrap_title <- function(x, width = 22) {
  paste(strwrap(x, width = width), collapse = "\n")
}

truncate_title <- function(x, max_chars = 20) {
  if (nchar(x) > max_chars) {
    paste0(substr(x, 1, max_chars), "...")
  } else {
    x
  }
}

plot_gene_heatmap_by_group <- function(df, gene_name, gbm_datasets, scale_rows = "none",
                                       cluster_rows = TRUE, cluster_cols = FALSE,
                                       na_col = "grey85", display_numbers = FALSE,
                                       cellwidth = NA, cellheight = NA, silent = FALSE) {
  
  sub_df <- df %>%
    filter(gene == gene_name) %>%
    select(-gene)
  
  mat <- sub_df %>%
    select(-file_name) %>%
    column_to_rownames("dataset") %>%
    as.matrix()
  
  if (nrow(mat) == 0) {
    warning("No data found for gene: ", gene_name)
    return(NULL)
  }
  
  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    message(n_na, " NA value(s) found — will be greyed out, not imputed")
  }
  
  is_gbm <- sapply(rownames(mat), function(fn) {
    any(sapply(gbm_datasets, function(ds) grepl(ds, fn, fixed = TRUE)))
  })
  
  annotation_row <- data.frame(Group = ifelse(is_gbm, "GBM", "non-GBM"))
  rownames(annotation_row) <- rownames(mat)
  
  annotation_colors <- list(
    Group = c(GBM = "#003366", `non-GBM` = "#D4A017")
  )
  
  row_order <- NULL
  if (isTRUE(cluster_rows)) {
    mat_for_clustering <- mat
    mat_for_clustering[is.na(mat_for_clustering)] <- 0
    
    gbm_names <- rownames(mat)[is_gbm]
    nongbm_names <- rownames(mat)[!is_gbm]
    
    order_within <- function(names_subset) {
      if (length(names_subset) <= 1) return(names_subset)
      sub_mat <- mat_for_clustering[names_subset, , drop = FALSE]
      d <- dist(sub_mat)
      hc <- hclust(d)
      names_subset[hc$order]
    }
    
    gbm_ordered <- order_within(gbm_names)
    nongbm_ordered <- order_within(nongbm_names)
    
    row_order <- c(gbm_ordered, nongbm_ordered)
    mat <- mat[row_order, , drop = FALSE]
    annotation_row <- annotation_row[row_order, , drop = FALSE]
  }
  
  pheatmap(
    mat,
    scale = scale_rows,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    na_col = na_col,
    display_numbers = display_numbers,
    number_format = "%.1f",
    angle_col = 90,
    border_color = NA,
    fontsize_col = 9,
    fontsize_row = 8,
    fontsize = 12,
    cellwidth = cellwidth,
    cellheight = cellheight,
    silent = silent,
    main = wrap_title(gene_name)
  )
}


plot_gene_heatmap_by_type <- function(df, gene_name, gbm_datasets, scale_rows = "none",
                                      cluster_rows = TRUE, cluster_cols = FALSE,
                                      na_col = "grey85", display_numbers = FALSE,
                                      cellwidth = NA, cellheight = NA, silent = FALSE) {
  
  sub_df <- df %>%
    filter(gene == gene_name) %>%
    select(-gene)
  
  mat <- sub_df %>%
    select(-file_name) %>%
    column_to_rownames("dataset") %>%
    as.matrix()
  
  if (nrow(mat) == 0) {
    warning("No data found for gene: ", gene_name)
    return(NULL)
  }
  
  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    message(n_na, " NA value(s) found — will be greyed out, not imputed")
  }
  
  is_gbm <- sapply(rownames(mat), function(fn) {
    any(sapply(gbm_datasets, function(ds) grepl(ds, fn, fixed = TRUE)))
  })
  
  annotation_row <- data.frame(Group = ifelse(is_gbm, "GBM", "non-GBM"))
  rownames(annotation_row) <- rownames(mat)
  
  annotation_colors <- list(
    Group = c(GBM = "#003366", `non-GBM` = "#D4A017")
  )
  
  row_order <- NULL
  if (isTRUE(cluster_rows)) {
    mat_for_clustering <- mat
    mat_for_clustering[is.na(mat_for_clustering)] <- 0
    
    gbm_names <- rownames(mat)[is_gbm]
    nongbm_names <- rownames(mat)[!is_gbm]
    
    order_within <- function(names_subset) {
      if (length(names_subset) <= 1) return(names_subset)
      sub_mat <- mat_for_clustering[names_subset, , drop = FALSE]
      d <- dist(sub_mat)
      hc <- hclust(d)
      names_subset[hc$order]
    }
    
    gbm_ordered <- order_within(gbm_names)
    nongbm_ordered <- order_within(nongbm_names)
    
    row_order <- c(gbm_ordered, nongbm_ordered)
    mat <- mat[row_order, , drop = FALSE]
    annotation_row <- annotation_row[row_order, , drop = FALSE]
  }
  
  pheatmap(
    mat,
    scale = scale_rows,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    na_col = na_col,
    display_numbers = display_numbers,
    number_format = "%.1f",
    angle_col = 90,
    border_color = NA,
    fontsize_col = 9,
    fontsize_row = 8,
    fontsize = 12,
    cellwidth = cellwidth,
    cellheight = cellheight,
    silent = silent,
    main = wrap_title(gene_name)
  )
}

plot_drug_heatmap_update <- function(df, drug_name, gbm_datasets, scale_rows = "none",
                                     cluster_rows = TRUE, cluster_cols = FALSE,
                                     na_col = "grey85", display_numbers = FALSE,
                                     cellwidth = NA, cellheight = NA, silent = FALSE) {
  
  # Filter for the specific drug
  sub_df <- df %>%
    filter(drug == drug_name) %>%
    select(-drug)
  
  mat <- sub_df %>%
    select(-file_name) %>%        # remove file_name if you don't want it in the heatmap
    column_to_rownames("dataset") %>%
    as.matrix()
  
  if (nrow(mat) == 0) {
    warning("No data found for drug: ", drug_name)
    return(NULL)
  }
  
  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    message(n_na, " NA value(s) found — will be greyed out, not imputed")
  }
  
  # Flag GBM by checking if file_name contains a GBM dataset ID
  is_gbm <- sapply(rownames(mat), function(fn) {
    any(sapply(gbm_datasets, function(ds) grepl(ds, fn, fixed = TRUE)))
  })
  
  annotation_row <- data.frame(Group = ifelse(is_gbm, "GBM", "non-GBM"))
  rownames(annotation_row) <- rownames(mat)
  
  annotation_colors <- list(
    Group = c(GBM = "#003366", `non-GBM` = "#D4A017")
  )
  
  # --- Cluster within each group using NA-safe matrix, then apply order to real mat ---
  row_order <- NULL
  if (isTRUE(cluster_rows)) {
    mat_for_clustering <- mat
    mat_for_clustering[is.na(mat_for_clustering)] <- 0
    
    gbm_names <- rownames(mat)[is_gbm]
    nongbm_names <- rownames(mat)[!is_gbm]
    
    order_within <- function(names_subset) {
      if (length(names_subset) <= 1) return(names_subset)
      sub_mat <- mat_for_clustering[names_subset, , drop = FALSE]
      d <- dist(sub_mat)
      hc <- hclust(d)
      names_subset[hc$order]
    }
    
    gbm_ordered <- order_within(gbm_names)
    nongbm_ordered <- order_within(nongbm_names)
    
    row_order <- c(gbm_ordered, nongbm_ordered)
    mat <- mat[row_order, , drop = FALSE]
    annotation_row <- annotation_row[row_order, , drop = FALSE]
  }
  pheatmap(
    mat,
    scale = scale_rows,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    na_col = na_col,
    display_numbers = display_numbers,
    number_format = "%.1f",
    angle_col = 90,        # added
    border_color = NA,     # added
    fontsize_col = 9,       # was 14
    fontsize_row = 8,       # added
    fontsize = 12,
    cellwidth = cellwidth,
    cellheight = cellheight,
    silent = silent,
    main = truncate_title(drug_name)
  )
}


plot_drug_heatmap_type <- function(df, drug_name, gbm_datasets, scale_rows = "none",
                                   cluster_rows = TRUE, cluster_cols = FALSE,
                                   na_col = "grey85", display_numbers = FALSE,
                                   cellwidth = NA, cellheight = NA, silent = FALSE) {
  
  sub_df <- df %>%
    filter(drug == drug_name) %>%
    select(-drug)
  
  mat <- sub_df %>%
    select(-file_name) %>%
    column_to_rownames("dataset") %>%
    as.matrix()
  
  if (nrow(mat) == 0) {
    warning("No data found for drug: ", drug_name)
    return(NULL)
  }
  
  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    message(n_na, " NA value(s) found — will be greyed out, not imputed")
  }
  
  is_gbm <- sapply(rownames(mat), function(fn) {
    any(sapply(gbm_datasets, function(ds) grepl(ds, fn, fixed = TRUE)))
  })
  
  annotation_row <- data.frame(Group = ifelse(is_gbm, "GBM", "non-GBM"))
  rownames(annotation_row) <- rownames(mat)
  
  annotation_colors <- list(
    Group = c(GBM = "#003366", `non-GBM` = "#D4A017")
  )
  
  row_order <- NULL
  if (isTRUE(cluster_rows)) {
    mat_for_clustering <- mat
    mat_for_clustering[is.na(mat_for_clustering)] <- 0
    
    gbm_names <- rownames(mat)[is_gbm]
    nongbm_names <- rownames(mat)[!is_gbm]
    
    order_within <- function(names_subset) {
      if (length(names_subset) <= 1) return(names_subset)
      sub_mat <- mat_for_clustering[names_subset, , drop = FALSE]
      d <- dist(sub_mat)
      hc <- hclust(d)
      names_subset[hc$order]
    }
    
    gbm_ordered <- order_within(gbm_names)
    nongbm_ordered <- order_within(nongbm_names)
    
    row_order <- c(gbm_ordered, nongbm_ordered)
    mat <- mat[row_order, , drop = FALSE]
    annotation_row <- annotation_row[row_order, , drop = FALSE]
  }
  pheatmap(
    mat,
    scale = scale_rows,
    cluster_rows = FALSE,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    na_col = na_col,
    display_numbers = display_numbers,
    number_format = "%.1f",
    angle_col = 90,
    border_color = NA,
    fontsize_col = 9,
    fontsize_row = 8,        # was 12 — now matches the other three
    fontsize = 12,
    cellwidth = cellwidth,
    cellheight = cellheight,
    silent = silent,
    main = truncate_title(drug_name)
  )
}


# Draw a heatmap at its natural size (fixed cell size) so cells and labels
# keep their proportions, instead of stretching a fixed-size PNG into a box.
render_heatmap_image <- function(plot_fn, ..., cellwidth, cellheight = 14, alt = "Heatmap", res = 192) {
  outfile <- tempfile(fileext = ".png")
  
  ph <- tryCatch(
    suppressMessages(suppressWarnings(
      plot_fn(..., cellwidth = cellwidth, cellheight = cellheight, silent = TRUE)
    )),
    error = function(e) NULL
  )
  
  if (is.null(ph)) {
    grDevices::png(outfile, width = 6, height = 1.5, units = "in", res = res)
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No single-cell data for this selection", col = "grey40")
    grDevices::dev.off()
    return(list(src = outfile, contentType = "image/png", alt = alt,
                style = "width: 576px; max-width: 100%; height: auto;"))
  }
  
  # Measure the finished heatmap (inches) on a null device
  grDevices::pdf(NULL)
  dev_id <- grDevices::dev.cur()
  size <- tryCatch(
    c(grid::convertWidth(sum(ph$gtable$widths), "in", valueOnly = TRUE),
      grid::convertHeight(sum(ph$gtable$heights), "in", valueOnly = TRUE)),
    error = function(e) c(NA_real_, NA_real_)
  )
  grDevices::dev.off(dev_id)
  if (anyNA(size) || any(size < 1)) size <- c(8, 7)
  size <- pmin(size + 0.3, c(18, 18))
  
  grDevices::png(outfile, width = size[1], height = size[2], units = "in", res = res)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  grDevices::dev.off()
  
  # Show at 96 px per inch (sharp on high-DPI screens), shrinking to fit the card
  list(
    src = outfile,
    contentType = "image/png",
    alt = alt,
    style = sprintf("width: %dpx; max-width: 100%%; height: auto; display: block; margin: 0 auto;",
                    round(size[1] * 96))
  )
}

get_cutpoint <- function(final_data, variable, method = "Median Cut") {
  
  if (method == "Median Cut") {
    return(median(final_data[[variable]], na.rm = TRUE))
  }
  
  if (method == "Maxstat Cut-off") {
    fmla <- as.formula(paste0("Surv(os_time, os_status) ~ ", variable))
    
    mstat <- maxstat_test(
      fmla,
      data = final_data,
      smethod = "LogRank",
      pmethod = "HL",
      abseps = 0.05
    )
    
    return(as.numeric(mstat@estimates$estimate$cutpoint))
  }
}

run_survival_plot_single_drug <- function(seq_final, Clinical_Data, drug_name, cohort_name, cutoff_method) {
  
  target_gene_ls <- drug_gene_ls[[drug_name]]
  
  filter_data <- tibble::add_column(
    as.data.frame(data_PC(seq_final, target_gene_ls, drug_name)$filter_data[, "PC_1"]),
    patient_barcode = seq_final[, "patient_barcode"]
  )
  colnames(filter_data)[1] <- "PC_1"
  
  final_data <- left_join(filter_data, Clinical_Data, by = "patient_barcode")
  
  cutpoint_1 <- get_cutpoint(final_data, variable = "PC_1", method = cutoff_method)
  
  df_sar_group <- final_data %>%
    dplyr::select(os_time, os_status, PC_1) %>%
    dplyr::mutate(
      drug_1_group = if_else(PC_1 <= cutpoint_1, "Low", "High"),
      drug_1_group = factor(drug_1_group, levels = c("Low", "High"))
    )
  
  surv.model <- survfit(Surv(os_time, os_status) ~ drug_1_group, data = df_sar_group)
  
  p <- ggsurvplot(
    surv.model,
    risk.table = TRUE,
    data = df_sar_group,
    pval = TRUE,
    ggtheme = nice_surv_theme,
    risk.table.fontsize = 4,
    risk.table.height = 0.25,
    xlab = "Overall survival (months)",
    break.time.by = 12,
    xlim = c(0, 60),
    legend.title = "Group",
    legend.labs = c("Low", "High")
  )
  
  return(p)
}

data_PC = function(data, target_gene_ls,drug_list){
  #Filter the gene expression data
  filter_data = data %>% 
    select(any_of(unlist(target_gene_ls,use.names = FALSE)))
  
  loading = list()
  
  for (i in seq_along(drug_list)){
    pca = prcomp(filter_data, center = TRUE, scale = TRUE)
    if (mean(pca$rotation[, 1] > 0) < 0.5)
    {
      pca$rotation <- pca$rotation * (-1)
      pca$x <- pca$x * (-1)
    }
    pca_value = unname(pca$x[,1])
    filter_data[[paste0("PC_", i)]] = pca_value
    loading[[drug_list[i]]] = as.data.frame(pca$rotation)
  }
  result = list(filter_data = filter_data,loading = loading)
  return(result)
}




ui <- page_fluid(
  theme = bs_theme(
    bootswatch = "flatly",
    base_font = font_google("Inter"),
    version = 5
  ),
  
  tags$head(
    tags$script(HTML("
      function tutorialZoom(img) {
        var box = document.getElementById('tutorial-zoom');
        box.querySelector('img').src = img.src;
        box.style.display = 'flex';
      }
      document.addEventListener('keydown', function(e) {
        if (e.key === 'Escape') {
          var box = document.getElementById('tutorial-zoom');
          if (box) box.style.display = 'none';
        }
      });
    ")),
    tags$style(HTML("
      .app-title {
        text-align: center;
        font-size: 30px;
        font-weight: 700;
        margin-top: 12px;
        margin-bottom: 4px;
      }
      .app-subtitle {
        text-align: center;
        color: #666;
        font-size: 15px;
        margin-bottom: 18px;
      }
      .card {
        border-radius: 16px !important;
        box-shadow: 0 4px 12px rgba(0,0,0,0.08);
        border: none;
      }
      th.has-tip {
        cursor: help;
        text-decoration: underline dotted;
        text-underline-offset: 3px;
      }
      th.has-tip::after {
        content: ' \\24D8';
        color: #18bc9c;
        font-weight: 400;
      }
      .tooltip-inner {
        max-width: 340px;
        text-align: left;
        font-size: 13px;
        line-height: 1.4;
      }
      .select-demo {
        display: flex; flex-wrap: wrap; gap: 16px; margin: 8px 0 14px 0;
      }
      .select-demo-step { flex: 1 1 220px; }
      .select-demo-label { font-weight: 600; margin-bottom: 6px; }
      .select-demo-box {
        border: 1px solid #ced4da; border-radius: 6px; padding: 6px 10px;
        background: #fff; display: flex; align-items: center; min-height: 36px;
        box-shadow: 0 0 0 3px rgba(24,188,156,0.15);
      }
      .select-demo-box .ghost { color: #aaa; }
      .select-demo-box .cursor { color: #2c3e50; animation: blink 1s step-end infinite; }
      .select-demo-box .caret { margin-left: auto; color: #888; }
      .select-demo-list {
        border: 1px solid #ced4da; border-top: none; border-radius: 0 0 6px 6px;
        background: #fff; font-size: 14px;
      }
      .select-demo-list div { padding: 4px 10px; }
      .select-demo-list .hit { background: #e8f8f4; }
      .select-demo-cap { font-size: 13px; color: #666; margin-top: 6px; }
      @keyframes blink { 50% { opacity: 0; } }
      .tutorial-step {
        display: flex; gap: 12px; align-items: flex-start; margin-bottom: 12px;
      }
      .tutorial-step .num {
        flex: 0 0 30px; height: 30px; border-radius: 50%;
        background: #2c3e50; color: #fff; font-weight: 700;
        display: flex; align-items: center; justify-content: center;
      }
      .tutorial-note {
        background: #f4f6f8; border-left: 4px solid #18bc9c;
        padding: 10px 14px; border-radius: 6px; font-size: 14px;
      }
      .featured-btns .btn { margin: 0 6px 6px 0; }
      .intro-lead { font-size: 16px; line-height: 1.55; max-width: 1000px; }
      .intro-stats { display: flex; flex-wrap: wrap; gap: 12px; margin: 14px 0 20px 0; }
      .intro-stat {
        flex: 1 1 150px; background: #f4f6f8; border-radius: 10px; padding: 12px 14px;
        border-top: 3px solid #18bc9c;
      }
      .intro-stat .val { font-size: 24px; font-weight: 700; color: #2c3e50; line-height: 1.1; }
      .intro-stat .lbl { font-size: 13px; color: #666; margin-top: 2px; }
      .intro-evidence { display: flex; flex-wrap: wrap; gap: 12px; }
      .intro-ev {
        flex: 1 1 260px; border: 1px solid #e3e6ea; border-radius: 10px; padding: 12px 14px;
        font-size: 14px; color: #444;
      }
      .intro-ev .ev-title { font-weight: 700; color: #2c3e50; margin-bottom: 4px; }
      .intro-ev .ev-title .fa, .intro-ev .ev-title svg { color: #18bc9c; }
      .intro-tabs { display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 10px 22px; }
      .intro-tab { font-size: 14px; color: #555; }
      .intro-tab a { font-weight: 700; font-size: 15px; }
      .intro-cta { margin-top: 18px; }
      .intro-cta .btn { margin-right: 8px; }
      .tab-help { font-size: 13px; color: #666; margin: 4px 0 0 0; }
      .intro-credit { font-size: 13px; color: #666; background: #f8f9fa; text-align: center; }
      .intro-credit a { color: #2c3e50; text-decoration: none; border-bottom: 1px dotted #2c3e50; }
      .intro-credit a:hover { color: #18bc9c; border-bottom-color: #18bc9c; }
      .drug-summary { font-size: 14px; padding-top: 2px; }
      .drug-summary .lbl { font-size: 12px; color: #888; text-transform: uppercase; letter-spacing: 0.03em; }
      .tutorial-heading {
        font-weight: 700; color: #2c3e50; margin: 6px 0 10px 0;
        display: flex; align-items: center; gap: 8px;
      }
      .tutorial-heading .fa, .tutorial-heading svg { color: #18bc9c; }
      .tutorial-bullets { margin: 6px 0 4px 0; padding-left: 18px; }
      .tutorial-bullets li { margin-bottom: 3px; }
      .evidence-mini td, .evidence-mini th { font-size: 13px; vertical-align: top; }
      .evidence-mini td:first-child { color: #666; white-space: nowrap; }
      #topCombTable tbody tr { cursor: pointer; }
      .tutorial-step > div:last-child { flex: 1; min-width: 0; }
      .tutorial-shot { margin: 10px 0 6px 0; }
      .tutorial-shot img {
        max-width: 100%; max-height: 420px; border: 1px solid #dee2e6;
        border-radius: 8px; box-shadow: 0 2px 6px rgba(0,0,0,0.08);
      }
      .tutorial-shot figcaption { font-size: 13px; color: #666; margin-top: 4px; }
      .tutorial-shot img { cursor: zoom-in; }
      .tutorial-gallery { display: flex; flex-wrap: wrap; gap: 14px; }
      .tutorial-gallery .tutorial-shot { flex: 1 1 280px; min-width: 0; }
      .tutorial-gallery .tutorial-shot img { width: 100%; max-height: 260px; object-fit: contain; background: #fff; }
      #tutorial-zoom {
        display: none; position: fixed; inset: 0; z-index: 2000;
        background: rgba(0,0,0,0.75); cursor: zoom-out;
        align-items: center; justify-content: center; padding: 30px;
      }
      #tutorial-zoom img { max-width: 95vw; max-height: 92vh; border-radius: 8px; background: #fff; }
      .tutorial-shot-missing {
        border: 2px dashed #ced4da; border-radius: 8px; background: #f8f9fa;
        color: #888; padding: 22px; text-align: center; font-size: 14px;
      }
    "))
  ),
  
  div(id = "tutorial-zoom", onclick = "this.style.display='none'", tags$img(alt = "Enlarged screenshot")),
  
  div(class = "app-title", "TACTIC-GBM: Tumor microenvironment – Aware Combination Therapy Index and Characterization for glioblastoma"),
  div(class = "app-subtitle", "Interactive exploration of genes, single drugs, doublets, and single-cell patterns"),
  
  navset_card_tab(
    id = "main_tabs",
    full_screen = TRUE,
    
    nav_panel(
      "Introduction",
      card(
        card_header("Welcome"),
        card_body(
          p(
            "TACTIC-GBM helps you explore candidate genes, therapies, and drug combinations",
            "for glioblastoma across TCGA and CGGA cohorts, and links those findings to",
            "single-cell expression patterns and drug-gene target databases."
          ),
          
          tags$hr(),
          
          tutorial_heading("list-check", "What you can do here"),
          
          tags$ul(
            tags$li(
              strong("Single Gene Analysis — "),
              "Select any candidate gene and view Kaplan-Meier survival curves stratified by",
              "expression level (Low/High) across TCGA and CGGA325 cohorts."
            ),
            tags$li(
              strong("Single Drug Analysis — "),
              "Explore PCA-based survival stratification for individual drugs, using their",
              "target-gene signatures. View top gene loadings for each cohort alongside",
              "survival curves, with a choice of Median or Maxstat cut-off methods."
            ),
            tags$li(
              strong("Top Combinations — "),
              "Browse a ranked table of candidate drug-pair (doublet) combinations, scored by",
              "clinical association, tumor microenvironment relevance, literature support, and",
              "combined rank-aggregation metrics. Filter, sort, and export results directly."
            ),
            tags$li(
              strong("Doublet Survival Analysis — "),
              "Pick a seed drug and a partner drug to view each drug's PC1-based (Maxstat) Low/High",
              "survival stratification side by side in TCGA and CGGA325."
            ),
            tags$li(
              strong("Single Cell Heatmap — "),
              "Visualize percent-expressing heatmaps by cell group and cell type, for either a",
              "selected drug's target genes or a single gene of interest, with GBM vs. non-GBM",
              "dataset annotations."
            ),
            tags$li(
              strong("Drug-Gene — "),
              "Look up any drug to see its annotated target genes across four source databases",
              "(DGIdb, DrugBank, DsigDB, DREIMT), and view Euler diagrams showing target-gene",
              "overlap across every combination of three sources."
            )
          ),
          
          tags$hr(),
          
          tutorial_heading("layer-group", "Three kinds of evidence"),
          div(
            class = "intro-evidence",
            div(class = "intro-ev",
                div(class = "ev-title", icon("heart-pulse"), " Clinical association"),
                "Each drug's target genes are summarized as one signature (PC1). A Cox model tests whether the two",
                " drugs' signatures are linked to overall survival in TCGA and CGGA325, and the two cohorts are combined",
                " with Fisher's method."),
            div(class = "intro-ev",
                div(class = "ev-title", icon("circle-nodes"), " Tumor microenvironment"),
                "Single-cell data show which cell populations (malignant, immune, vascular, neuronal, stromal)",
                " express each drug's targets, and how specifically."),
            div(class = "intro-ev",
                div(class = "ev-title", icon("book-open"), " Literature"),
                "PubMed abstracts are scored for relevance to the drug, its targets and glioblastoma.")
          ),
          p(style = "font-size: 14px; color: #555; margin-top: 6px;",
            "The three rankings are merged with Robust Rank Aggregation into the TACTIC score (RRA P-value):",
            " lower means a doublet ranks near the top on all three at once."),
          
          tags$hr(),
          
          tutorial_heading("compass", "Explore the app"),
          div(
            class = "intro-tabs",
            div(class = "intro-tab", actionLink("intro_go_top", "Top Combinations"),
                div("The 81 ranked doublets with all scores. Search, filter, export, and download an evidence package or PubMed IDs.")),
            div(class = "intro-tab", actionLink("intro_go_doublet", "Doublet Survival Analysis"),
                div("Survival curves for both drugs of a doublet, side by side in TCGA and CGGA325.")),
            div(class = "intro-tab", actionLink("intro_go_single_drug", "Single Drug Analysis"),
                div("One drug's target-gene signature: survival curves and the genes that drive it (PC1 loadings).")),
            div(class = "intro-tab", actionLink("intro_go_single_gene", "Single Gene Analysis"),
                div("Survival by High/Low expression of any single target gene.")),
            div(class = "intro-tab", actionLink("intro_go_heatmap", "Single Cell Heatmap"),
                div("Which cell populations express a drug's targets, or a single gene, across GBM and non-GBM datasets.")),
            div(class = "intro-tab", actionLink("intro_go_drug_gene", "Drug-Gene"),
                div("A drug's target genes in DGIdb, DrugBank, DSigDB and DREIMT."))
          ),
          
          div(
            class = "intro-cta",
            actionButton("intro_btn_tutorial", "Start the tutorial", icon = icon("graduation-cap"), class = "btn-primary"),
            actionButton("intro_btn_top", "See the top doublets", icon = icon("ranking-star"), class = "btn-outline-primary")
          ),
          
          tags$hr(),
          
          p(
            style = "color: #888; font-size: 13px;",
            "Use the tabs above to navigate between analyses. Most tabs include selectable",
            "inputs (genes, drugs, cut-off methods) that update the plots and tables in real time.",
            "New here? Start with the ", actionLink("intro_to_tutorial", "Tutorial"), " tab."
          )
        ),
        card_footer(
          class = "intro-credit",
          icon("flask"),
          " Created by the ",
          tags$a(href = "https://medicine.utah.edu/internal-medicine/epidemiology/chen-lab",
                 target = "_blank", rel = "noopener", strong("Chen Lab")),
          ", Division of Epidemiology, Department of Internal Medicine, University of Utah School of Medicine."
        )
      )
    ),
    
    nav_panel(
      "Tutorial",
      card(
        card_header("Tutorial"),
        card_body(
          tutorial_heading("route", "How to use TACTIC-GBM"),
          p("TACTIC-GBM ranks drug doublets for glioblastoma using three independent kinds of evidence,",
            "then lets you inspect each one. A typical session follows four steps:"),
          div(class = "tutorial-note", style = "margin-bottom: 14px;",
              icon("circle-info"), " Want to know how to select a drug or gene? See ",
              tags$a(href = "#tut-choose",
                     onclick = "document.getElementById('tut-choose').scrollIntoView({behavior: 'smooth', block: 'start'}); return false;",
                     strong("How to choose a drug or gene")),
              " below."),
          div(class = "tutorial-step", div(class = "num", "1"), div(
            strong("Find candidates. "),
            "Open ", actionLink("tut_go_top", "Top Combinations"), " to browse the 81 prioritized doublets,",
            "sorted by RRA P-value (lower is better). Hover any column name for its definition.",
            tags$ul(
              class = "tutorial-bullets",
              tags$li("To look at one drug, type its name in the ", strong("Search"), " box at the top right of",
                      " the table (for example ", tags$code("Sag1.3"), "). Only doublets containing that drug stay listed."),
              tags$li("To filter a single column, type in the box directly under that column's name."),
              tags$li("At the bottom of the tab you can download an evidence package for a doublet,",
                      " or the PubMed IDs for a drug (see ", em("How to download evidence and PMIDs"), " below).")
            ),
            tutorial_gallery(
              c("step1_top_combination.png", "step1_evidence.png"),
              c("Top Combinations filtered to Sag1.3 with the Search box at the top right; hover a column name for its explanation.",
                "Bottom of the tab: download the evidence package for one doublet (left) or PubMed IDs for any drugs (right).")
            )
          )),
          div(class = "tutorial-step", div(class = "num", "2"), div(
            strong("Check clinical evidence. "),
            "In ", actionLink("tut_go_doublet", "Doublet Survival Analysis"), ", select the two drugs to see",
            "how each drug's target-gene signature (PC1) stratifies overall survival in TCGA (top) and CGGA325 (bottom).",
            "Use ", actionLink("tut_go_single_drug", "Single Drug Analysis"), " for one drug's PC1 loadings,",
            "and ", actionLink("tut_go_single_gene", "Single Gene Analysis"), " for an individual target gene.",
            tutorial_gallery(
              c("step2_doublet.png", "step2_singledrug.png", "step2_singlegene.png"),
              c("Doublet Survival Analysis: GW-559090 + CHEMBL1321401, TCGA on top and CGGA325 below.",
                "Single Drug Analysis: type a name (here PD-184352) to see its PC1 loadings.",
                "Single Gene Analysis: survival by High/Low expression of one gene.")
            )
          )),
          div(class = "tutorial-step", div(class = "num", "3"), div(
            strong("Check the tumor microenvironment. "),
            "In ", actionLink("tut_go_heatmap", "Single Cell Heatmap"), ", see which cell populations express",
            "the drug's target genes, by merged cell group (malignant, immune, vascular, neuronal, stromal)",
            "and by the original studies' cell types. GBM datasets are labelled in navy, non-GBM in gold.",
            tutorial_gallery(
              c("step3_sc_bydrug.png", "step3_sc_bygene.png"),
              c("By drug: percentage of cells expressing GW-559090's target genes, by cell group (left) and cell type (right).",
                "By gene: the same view for a single gene, here ITGA5.")
            )
          )),
          div(class = "tutorial-step", div(class = "num", "4"), div(
            strong("Check the targets. "),
            "In ", actionLink("tut_go_drug_gene", "Drug-Gene"), ", look up which databases (DGIdb, DrugBank,",
            "DSigDB, DREIMT) annotate the drug's target genes."
          )),
          
          tags$hr(),
          
          div(id = "tut-choose", style = "scroll-margin-top: 20px;"),
          tutorial_heading("hand-pointer", "How to choose a drug or gene"),
          p("Every Drug and Gene box is a searchable list with thousands of names, so typing is faster than scrolling."),
          div(
            class = "select-demo",
            div(class = "select-demo-step",
                div(class = "select-demo-label", "1. Click the box"),
                div(class = "select-demo-box", span(class = "ghost", "Gw-559090"), span(class = "caret", "▾")),
                div(class = "select-demo-cap", "The current name turns grey and the box is ready for typing.")),
            div(class = "select-demo-step",
                div(class = "select-demo-label", "2. Type part of the name"),
                div(class = "select-demo-box", span("rwj"), span(class = "cursor", "|"), span(class = "caret", "▾")),
                div(class = "select-demo-list", div(class = "hit", "Rwj-67657")),
                div(class = "select-demo-cap", "Matches appear as you type. Upper/lower case doesn't matter.")),
            div(class = "select-demo-step",
                div(class = "select-demo-label", "3. Pick it"),
                div(class = "select-demo-box", span("Rwj-67657"), span(class = "caret", "▾")),
                div(class = "select-demo-cap", "Click the name or press Enter. Plots and tables update automatically."))
          ),
          div(class = "tutorial-note",
              strong("Tips. "),
              tags$ul(
                class = "tutorial-bullets",
                tags$li("If the old name doesn't clear when you click, press ", tags$kbd("Backspace"),
                        " to delete it, then type."),
                tags$li("If you click away without picking anything, the previous selection comes back."),
                tags$li("No match? The drug or gene isn't in that tab's data. Try a shorter part of the name ",
                        "(for example ", tags$code("chembl1321401"), " instead of the full ", tags$code("Chembl:Chembl1321401"), ")."),
                tags$li("The drugs discussed in the manuscript are listed first, under ",
                        tags$em("Discussed in manuscript"), ".")
              )),
          
          tags$hr(),
          
          tutorial_heading("download", "How to download evidence and PMIDs"),
          p("Both downloads are at the bottom of the ", actionLink("tut_go_top_dl", "Top Combinations"), " tab."),
          tags$ul(
            class = "tutorial-bullets",
            tags$li(strong("Evidence package (one doublet). "),
                    "Click a row in the table, or pick the doublet in the ", em("Doublet"), " box, then click ",
                    strong("Download evidence package"), ". The CSV has one row per drug with all the doublet's scores,",
                    " each drug's target genes and source, its dominant cell group, and its PubMed IDs."),
            tags$li(strong("PubMed IDs (any drugs). "),
                    "Type one or more drug names in the ", em("Drugs"), " box to preview their literature, then click ",
                    strong("Download PMIDs"), ". The CSV has one row per article, with the drug's target genes,",
                    " its dominant cell group and a PubMed link.")
          ),
          
          tags$hr(),
          
          tutorial_heading("flask", "How to try a worked example"),
          p("GW-559090 (targets ITGA5, DrugBank) + CHEMBL1321401 (targets BMP2, DGIdb) is rank 1 in",
            "manuscript Table 1 (Clinical Association Score 3.81, RRA P-value < 0.001).",
            "Higher ITGA5 expression is associated with poorer survival and higher BMP2 with better survival,",
            "and both genes are expressed across malignant, vascular and myeloid cells."),
          p("Load this example into each tab:"),
          div(
            class = "featured-btns",
            actionButton("ex_doublet", "Survival curves", icon = icon("chart-line"), class = "btn-outline-primary btn-sm"),
            actionButton("ex_single_drug", "PC1 loadings (GW-559090)", icon = icon("table"), class = "btn-outline-primary btn-sm"),
            actionButton("ex_heatmap", "Single-cell heatmaps", icon = icon("table-cells"), class = "btn-outline-primary btn-sm"),
            actionButton("ex_gene", "Single gene survival", icon = icon("dna"), class = "btn-outline-primary btn-sm"),
            actionButton("ex_drug_gene", "Target genes", icon = icon("database"), class = "btn-outline-primary btn-sm")
          ),
          p(style = "font-size: 13px; color: #666;",
            "The other two doublets discussed in the manuscript, PD-184352 + SAG1.3 and RWJ-67657 + AS-1409,",
            "appear at the top of the drug lists and as quick-pick buttons in Doublet Survival Analysis."),
          
          div(class = "tutorial-note",
              strong("Cohorts. "),
              "Primary IDH-wildtype GBM only: 141 TCGA and 74 CGGA325 patients. Survival is overall survival in months.")
        )
      )
    ),
    
    nav_panel(
      "Drug-Gene",
      
      card(
        card_header("Drug–Gene Overview"),
        card_body(
          fluidRow(
            column(
              width = 5,
              selectizeInput(
                inputId = "drug_gene_selected",
                label = "Select a drug",
                choices = drug_gene_tab_choices,
                selected = drug_gene_tab_choices[1],
                options = click_to_type_options("Type a drug name...")
              ),
              DTOutput("drugGeneTable")
            ),
            
            column(
              width = 7,
              p(strong("Target-Gene Overlap Across Databases"), style = "margin-bottom: 4px;"),
              p("Each diagram compares target genes for the selected drug across three of the four sources.",
                style = "font-size: 13px; color: #666;"),
              
              fluidRow(
                column(6, plotOutput("drugGeneVenn1", height = "260px")),
                column(6, plotOutput("drugGeneVenn2", height = "260px"))
              ),
              fluidRow(
                column(6, plotOutput("drugGeneVenn3", height = "260px")),
                column(6, plotOutput("drugGeneVenn4", height = "260px"))
              )
            )
          )
        )
      )
    ),
    
    nav_panel(
      "Single Gene Analysis",
      layout_column_wrap(
        width = 1,
        card(
          card_header("Single Gene Survival"),
          card_body(
            p("Select a gene to visualize survival stratification across cohorts."),
            selectizeInput("gene", "Gene", choices = unique_genes,
                           options = click_to_type_options("Type a gene name...")),
            plotOutput("survivalPlot", height = "500px")
          )
        )
      )
    ),
    
    nav_panel(
      "Single Drug Analysis",
      layout_column_wrap(
        width = 1,
        fill = FALSE,
        heights_equal = "row",
        card(
          fill = FALSE,
          card_header("Single Drug Survival"),
          card_body(
            fill = FALSE,
            layout_columns(
              col_widths = c(5, 4, 3),
              fill = FALSE,
              selectizeInput("drug", "Drug", choices = single_drug_choices, width = "100%",
                             options = click_to_type_options("Type a drug name...")),
              radioButtons("cutoff_method", "Cut-off method",
                           choices = c("Median" = "Median Cut", "Maxstat (optimal)" = "Maxstat Cut-off"),
                           inline = TRUE),
              uiOutput("drug_summary")
            ),
            p(class = "tab-help",
              icon("circle-info"),
              " The drug's target genes are summarized as one signature score per patient (PC1).",
              " Patients are split into Low and High at the median, or at the Maxstat-optimal cut point.")
          )
        ),
        layout_columns(
          col_widths = c(6, 6),
          fill = FALSE,
          card(fill = FALSE, card_header("TCGA: PC1 loadings"), card_body(fill = FALSE, DTOutput("pcLoading_tcga"))),
          card(fill = FALSE, card_header("CGGA325: PC1 loadings"), card_body(fill = FALSE, DTOutput("pcLoading_cgga325")))
        ),
        p(class = "tab-help", style = "margin-top: -6px;",
          "Loadings show how much each target gene contributes to the signature. PC1 is oriented so most",
          " loadings are positive, so \"High\" means higher overall target-gene expression.",
          " Genes are sorted by the size of their PC1 loading."),
        card(
          fill = FALSE,
          full_screen = TRUE,
          card_header("Survival curves"),
          card_body(fill = FALSE, plotOutput("drugSurvivalPlot", height = "460px"))
        )
      )
    ),
    
    nav_panel(
      "Top Combinations",
      layout_column_wrap(
        width = 1,
        card(
          card_header("Candidate Doublets"),
          card_body(
            p("Ranked candidate combinations, sorted by RRA P-value (lower = stronger, more consistent support).",
              "Hover a column header for its definition. Type a drug name in the Search box (top right) to show only",
              "doublets containing it. Click a row to select it for the evidence package below."),
            DTOutput("topCombTable")
          )
        ),
        
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header(icon("file-zipper"), " Download evidence package (one doublet)"),
            card_body(
              selectizeInput("evidence_doublet", "Doublet", choices = doublet_pair_choices,
                             width = "100%",
                             options = click_to_type_options("Type a drug name...")),
              uiOutput("evidence_preview"),
              downloadButton("download_evidence", "Download evidence package (CSV)", class = "btn-primary")
            )
          ),
          card(
            card_header(icon("book-medical"), " Download PubMed IDs (any drugs)"),
            card_body(
              selectizeInput("pmid_drugs", "Drugs", choices = pmid_drug_choices, multiple = TRUE,
                             width = "100%",
                             options = list(placeholder = "Type one or more drug names...",
                                            plugins = list("remove_button"))),
              DTOutput("pmid_preview"),
              br(),
              downloadButton("download_pmids", "Download PMIDs (CSV)", class = "btn-primary")
            )
          )
        )
      )
    ),
    
    nav_panel(
      "Doublet Survival Analysis",
      layout_column_wrap(
        width = 1,
        card(
          card_header("Two-Drug Survival Comparison (Maxstat Cut-off)"),
          card_body(
            p("Each column shows one selected drug's Maxstat-based Low/High stratification;",
              "TCGA is shown on top, CGGA325 below."),
            fluidRow(
              column(6, selectizeInput("drug1_pc1", "Drug 1", choices = doublet_drug1_choices,
                                       options = click_to_type_options("Type a drug name..."))),
              column(6, selectizeInput("drug2_pc1", "Drug 2", choices = doublet_drug2_choices,
                                       options = click_to_type_options("Type a drug name...")))
            ),
            div(
              class = "featured-btns",
              span("Doublets discussed in the manuscript: ",
                   style = "font-size: 13px; color: #666; margin-right: 4px;"),
              lapply(seq_len(nrow(featured_doublets)), function(i) {
                actionButton(
                  paste0("featured_doublet_", i),
                  featured_doublets$label[i],
                  class = "btn-outline-secondary btn-sm",
                  title = sprintf("Table 1 rank %d, RRA P-value %s",
                                  featured_doublets$table1_rank[i], featured_doublets$rra_p[i])
                )
              })
            ),
            br(),
            plotOutput("pc1SurvivalPlot", height = "900px")
          )
        )
      )
    ),
    
    nav_panel(
      "Single Cell Heatmap",
      layout_column_wrap(
        width = 1,
        card(
          card_header("Mechanistic Insights — By Drug"),
          card_body(
            fluidRow(
              column(
                width = 6,
                selectizeInput("Drug1", "Drug", choices = heatmap_drug_choices,
                               options = click_to_type_options("Type a drug name..."))
              ),
              column(
                width = 6,
                div(
                  style = "margin-top: 26px; color: #666; font-size: 13px; line-height: 1.4;",
                  icon("circle-info", style = "margin-right: 4px;"),
                  "Each heatmap shows the percentage of cells expressing the selected drug's target genes."
                )
              )
            ),
            br(),
            fluidRow(
              column(6, card(card_header("Cell Group"), card_body(imageOutput("HeatmapPlot_by_group", height = "auto")))),
              column(6, card(card_header("Cell Type"), card_body(imageOutput("HeatmapPlot_by_type", height = "auto"))))
            ),
            uiOutput("drug1_full_name_note")
          )
        ),
        
        card(
          card_header("Mechanistic Insights — By Gene"),
          card_body(
            fluidRow(
              column(
                width = 6,
                selectizeInput("Gene1", "Gene", choices = unique(merge_by_group_gene$gene),
                               options = click_to_type_options("Type a gene name..."))
              ),
              column(
                width = 6,
                div(
                  style = "margin-top: 26px; color: #666; font-size: 13px; line-height: 1.4;",
                  icon("circle-info", style = "margin-right: 4px;"),
                  "Each heatmap shows the percentage of cells expressing the selected gene."
                )
              )
            ),
            br(),
            fluidRow(
              column(6, card(card_header("Cell Group"), card_body(imageOutput("GeneHeatmapPlot_by_group", height = "auto")))),
              column(6, card(card_header("Cell Type"), card_body(imageOutput("GeneHeatmapPlot_by_type", height = "auto"))))
            ),
            uiOutput("gene1_full_name_note")
          )
        )
      )
    ),
    nav_panel(
      "Contact",
      layout_column_wrap(
        width = 1,
        card(
          card_header("Send Us a Message"),
          card_body(
            p("Leave your name, email, and message below and our team will get back to you.",
              style = "color: #666; font-size: 13px;"),
            textInput("sender_name", "Your Name", placeholder = "Jane Doe"),
            textInput("sender_email", "Your Email", placeholder = "jane@example.com"),
            textInput("subject", "Subject", placeholder = "Question about the app"),
            textAreaInput("message", "Message", rows = 6, width = "100%",
                          placeholder = "Write your message here..."),
            actionButton("send_btn", "Send Message", class = "btn-primary",
                         icon = icon("paper-plane")),
            div(id = "contact_status", style = "margin-top: 10px; font-size: 14px;"),
            
            # Sends the form from the visitor's browser to Web3Forms
            tags$script(HTML(sprintf("
              document.addEventListener('click', async function (e) {
                const btn = e.target.closest('#send_btn');
                if (!btn) return;
                
                const val = id => (document.getElementById(id).value || '').trim();
                const status = document.getElementById('contact_status');
                const name = val('sender_name'), email = val('sender_email'),
                      subject = val('subject'), message = val('message');
                
                if (!email || !message) {
                  status.style.color = '#c0392b';
                  status.textContent = 'Please enter your email and a message.';
                  return;
                }
                if (!/^[^@\\s]+@[^@\\s]+\\.[^@\\s]+$/.test(email)) {
                  status.style.color = '#c0392b';
                  status.textContent = 'Please enter a valid email address.';
                  return;
                }
                
                btn.disabled = true;
                status.style.color = '#666';
                status.textContent = 'Sending...';
                
                try {
                  const res = await fetch('https://api.web3forms.com/submit', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
                    body: JSON.stringify({
                      access_key: '%s',
                      subject: 'TACTIC-GBM Contact Form: ' + (subject || '(no subject)'),
                      from_name: 'TACTIC-GBM App',
                      name: name || '(not provided)',
                      email: email,        // lets Ann click Reply to answer the sender
                      message: message
                    })
                  });
                  const out = await res.json();
                  if (out.success) {
                    status.style.color = '#18bc9c';
                    status.textContent = 'Message sent successfully! We will get back to you soon.';
                    ['sender_name', 'sender_email', 'subject', 'message'].forEach(function (id) {
                      const el = document.getElementById(id);
                      el.value = '';
                      el.dispatchEvent(new Event('change'));
                    });
                  } else {
                    status.style.color = '#c0392b';
                    status.textContent = 'Failed to send: ' + (out.message || 'unknown error');
                  }
                } catch (err) {
                  status.style.color = '#c0392b';
                  status.textContent = 'Failed to send. Please check your internet connection and try again.';
                } finally {
                  btn.disabled = false;
                }
              });
            ", web3forms_access_key)))
          )
        )
      )
    )
  )
)


####### Here is the part for Data Preparation ######


####### Here is the part for Shiny Server #######

server <- function(input, output, session) {
  
  geneData <- reactive({
    selected_gene <- input$gene
    
    list(tcga = gbm_tcga_seq_final_geneanlysis_primary,
         #cgga693 = CGGA_rna_seq_final_693_geneanlysis_primary,
         cgga325 = CGGA_rna_seq_final_325_geneanlysis_primary,
         gene = selected_gene)
  })
  
  
  ClinicalData <- reactive({
    selected_gene <- input$gene
    
    list(tcga = tcga_clinical_subtype,
         #cgga693 = cgga693_clinical_subtype,
         cgga325 = cgga325_clinical_subtype,
         gene = selected_gene)
  })
  
  
  selected_drug_gene_sets <- reactive({
    
    req(input$drug_gene_selected)
    
    selected_drug <- input$drug_gene_selected
    
    list(
      DGIdb = get_drug_genes_by_source(
        drug_name = selected_drug,
        source_name = "DGIdb",
        database_list = drug_db_ls_new
      ),
      
      DrugBank = get_drug_genes_by_source(
        drug_name = selected_drug,
        source_name = "DrugBank",
        database_list = drug_db_ls_new
      ),
      
      DsigDB = get_drug_genes_by_source(
        drug_name = selected_drug,
        source_name = "DsigDB",
        database_list = drug_db_ls_new
      ),
      
      Dreimt = get_drug_genes_by_source(
        drug_name = selected_drug,
        source_name = "DREIMT",
        database_list = drug_db_ls_new
      )
    )
  })
  
  
  drug_sets <- list(
    DGIdb    = names(drug_set_ls_grouped)[sapply(drug_set_ls_grouped, function(x) "DGIdb" %in% x)],
    DrugBank = names(drug_set_ls_grouped)[sapply(drug_set_ls_grouped, function(x) "DrugBank" %in% x)],
    DsigDB   = names(drug_set_ls_grouped)[sapply(drug_set_ls_grouped, function(x) "DsigDB" %in% x)],
    Dreimt   = names(drug_set_ls_grouped)[sapply(drug_set_ls_grouped, function(x) "DREIMT" %in% x)]
  )
  
  # Light color palette
  fill_colors <- c("#F6B8A2", "#8FD9D0", "#B8E0A0")
  
  # All combinations of 3 out of 4 sets
  combos <- combn(names(drug_sets), 3, simplify = FALSE)
  
  # Store plots in a list
  plot_list <- list()
  
  for (i in seq_along(combos)) {
    sets <- combos[[i]]
    s1 <- drug_sets[[sets[1]]]
    s2 <- drug_sets[[sets[2]]]
    s3 <- drug_sets[[sets[3]]]
    
    fit <- euler(c(
      setNames(length(s1), sets[1]),
      setNames(length(s2), sets[2]),
      setNames(length(s3), sets[3]),
      setNames(length(intersect(s1, s2)), paste(sets[1], sets[2], sep = "&")),
      setNames(length(intersect(s1, s3)), paste(sets[1], sets[3], sep = "&")),
      setNames(length(intersect(s2, s3)), paste(sets[2], sets[3], sep = "&")),
      setNames(length(Reduce(intersect, list(s1, s2, s3))), 
               paste(sets[1], sets[2], sets[3], sep = "&"))
    ), shape = "ellipse")
    
    plot_list[[i]] <- plot(fit,
                           fills = list(fill = fill_colors, alpha = 0.85),
                           labels = list(col = "grey20", font = 2, cex = 1),
                           quantities = list(type = "counts", cex = 0.9, col = "grey20", font = 2),
                           edges = list(col = "white", lwd = 2),
                           main = paste(sets, collapse = " vs ")
    )
  }
  
  output$survivalPlot <- renderPlot({
    
    plots <- list()
    
    for (cohort_name in c("tcga", "cgga325")) {
      
      
      gene_select <- geneData()$gene
      gene_Data <- geneData()[[cohort_name]]
      Clinical_Data <- ClinicalData()[[cohort_name]]
      
      gene_data_select <- gene_Data[, c("patient_barcode", gene_select)]
      final_data <- merge(gene_data_select, Clinical_Data, by = "patient_barcode")
      
      cutpoint_1 <- median(final_data[[gene_select]], na.rm = TRUE)
      
      df_sar_group <- final_data %>% 
        dplyr::mutate(gene_group = if_else(final_data[[gene_select]] <= cutpoint_1, "Low", "High")) 
      
      surv.model <- survfit(Surv(os_time, os_status) ~ gene_group, data = df_sar_group)
      
      p <- ggsurvplot(surv.model, pval = FALSE, risk.table = TRUE,
                      data = df_sar_group,
                      font.tickslab = c(18),
                      ggtheme = theme_classic2(base_size = 15), 
                      risk.table.fontsize = 10, 
                      tables.height = 0.5, 
                      xlab = "Survival time in months")
      
      p$plot <- p$plot + ggtitle(toupper(cohort_name))
      
      plots[[cohort_name]] <- p$plot
    }
    
    grid.arrange(plots$tcga, plots$cgga325, ncol = 2)
  })
  
  
  output$drugSurvivalPlot <- renderPlot({
    req(input$drug, input$cutoff_method)
    plots <- list()
    
    for (cohort_name in c("tcga", "cgga325")) {
      p <- run_survival_plot_single_drug(
        seq_final = geneData()[[cohort_name]],
        Clinical_Data = ClinicalData()[[cohort_name]],
        drug_name = input$drug,
        cohort_name = toupper(cohort_name),
        cutoff_method = input$cutoff_method
      )
      plots[[cohort_name]] <- p$plot + ggtitle(paste0(toupper(cohort_name), ": ", input$drug))
    }
    
    gridExtra::grid.arrange(grobs = plots, ncol = 2)
  })
  
  # PCA loadings for the selected drug, computed once per cohort
  drug_loadings <- reactive({
    req(input$drug)
    genes <- drug_gene_ls[[input$drug]]
    lapply(c(tcga = "tcga", cgga325 = "cgga325"), function(cohort) {
      res <- data_PC(geneData()[[cohort]], genes, input$drug)
      as.data.frame(res$loading[[input$drug]]) %>%
        tibble::rownames_to_column("Gene") %>%
        dplyr::select(Gene, dplyr::any_of(c("PC1", "PC2", "PC3"))) %>%
        dplyr::arrange(dplyr::desc(abs(PC1)))
    })
  })
  
  output$drug_summary <- renderUI({
    req(input$drug)
    genes <- sort(unique(unlist(drug_gene_ls[[input$drug]])))
    measured <- intersect(genes, colnames(geneData()$tcga))
    sources <- drug_set_ls_grouped[[input$drug]]
    div(
      class = "drug-summary",
      div(class = "lbl", "Target genes"),
      div(strong(length(genes)),
          if (length(measured) < length(genes)) span(sprintf(" (%d measured)", length(measured)), style = "color:#888;")),
      if (length(sources)) tagList(
        div(class = "lbl", style = "margin-top: 6px;", "Source"),
        div(lapply(sources, function(x) span(class = "badge text-bg-light", style = "margin-right: 4px;", x)))
      )
    )
  })
  
  render_loading_table <- function(tab) {
    n <- nrow(tab)
    dt <- datatable(
      tab,
      rownames = FALSE,
      options = list(
        dom = if (n > 8) "tp" else "t",
        pageLength = 8,
        ordering = FALSE,
        scrollX = TRUE,
        columnDefs = list(list(className = "dt-left", targets = 0))
      )
    ) %>%
      formatRound(columns = setdiff(names(tab), "Gene"), digits = 3)
    if ("PC1" %in% names(tab)) {
      dt <- dt %>%
        formatStyle("PC1", fontWeight = "bold",
                    color = styleInterval(0, c("#c0392b", "#1a7f37")))
    }
    dt
  }
  
  output$pcLoading_tcga <- DT::renderDT(render_loading_table(drug_loadings()$tcga))
  output$pcLoading_cgga325 <- DT::renderDT(render_loading_table(drug_loadings()$cgga325))
  
  output$topCombTable <- DT::renderDT({
    # Same row order as candidate_ranked, so a clicked row maps to its doublet
    tab <- candidate_ranked %>%
      dplyr::select(-rra_rank, -pair_label) %>%
      # Show P and -log10(P) side by side in place of the stored Fisher column
      dplyr::relocate(clinical_p, clinical_score, .after = Fisher_combined_tcga_cgga325) %>%
      dplyr::select(-Fisher_combined_tcga_cgga325) %>%
      dplyr::rename(dplyr::any_of(c(
        `Drug 1` = "Drug1",
        `Drug 1 Source` = "drug_set_drug1",
        `Drug 2` = "Drug2",
        `Drug 2 Source` = "drug_set_drug2",
        `Fisher-combined P-value` = "clinical_p",
        `Clinical Association Score` = "clinical_score",
        `TME Specificity Score` = "TME_score",
        `Literature Semantic Score` = "literature_score",
        `Reciprocal Rank Fusion Score` = "RRF_score",
        `RRA P-value (TACTIC Score)` = "RRA_score"
      )))
    
    # Header tooltips: Bootstrap tooltip when available, native title otherwise
    tooltip_js <- DT::JS(sprintf(
      "function(settings) {
         var tips = %s;
         var api = this.api();
         $(api.table().container()).find('thead th').each(function() {
           var t = $(this).text().trim();
           if (tips[t]) {
             $(this).attr('title', tips[t]).addClass('has-tip');
             if (window.bootstrap && bootstrap.Tooltip) {
               new bootstrap.Tooltip(this, {container: 'body', placement: 'top'});
             }
           }
         });
       }",
      jsonlite::toJSON(score_definitions, auto_unbox = TRUE)
    ))
    
    numeric_cols <- names(tab)[sapply(tab, is.numeric)]
    
    datatable(
      tab,
      rownames = FALSE,
      filter = "top",
      selection = "single",
      extensions = c("Buttons", "Scroller"),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel"),
        autoWidth = TRUE,
        initComplete = tooltip_js
      )
    ) %>%
      formatRound(
        columns = setdiff(numeric_cols, c("Fisher-combined P-value", "Clinical Association Score")),
        digits = 3
      ) %>%
      formatRound(columns = "Clinical Association Score", digits = 2) %>%
      formatSignif(columns = "Fisher-combined P-value", digits = 2)
  })
  
  ## Evidence package ---------------------------------------------------
  # Clicking a table row selects that doublet for download
  observeEvent(input$topCombTable_rows_selected, {
    i <- input$topCombTable_rows_selected
    if (length(i)) {
      updateSelectizeInput(session, "evidence_doublet", selected = candidate_ranked$pair_label[i[1]])
    }
  })
  
  selected_pair <- reactive({
    req(input$evidence_doublet)
    candidate_ranked[candidate_ranked$pair_label == input$evidence_doublet, , drop = FALSE][1, ]
  })
  
  evidence_table <- reactive(build_evidence_package(selected_pair()))
  
  output$evidence_preview <- renderUI({
    ev <- evidence_table()
    pr <- selected_pair()
    div(
      class = "evidence-preview",
      div(strong(ev$Doublet[1]), span(sprintf("  RRA rank %d of %d", pr$rra_rank, nrow(candidate_ranked)),
                                      style = "color: #888;")),
      tags$table(
        class = "table table-sm evidence-mini",
        tags$thead(tags$tr(tags$th(""), tags$th(ev$Drug[1]), tags$th(ev$Drug[2]))),
        tags$tbody(
          tags$tr(tags$td("Target genes"), tags$td(ev$`Target genes`[1]), tags$td(ev$`Target genes`[2])),
          tags$tr(tags$td("Dominant cell group"),
                  tags$td(ifelse(is.na(ev$`Dominant cell group`[1]), "n/a", ev$`Dominant cell group`[1])),
                  tags$td(ifelse(is.na(ev$`Dominant cell group`[2]), "n/a", ev$`Dominant cell group`[2]))),
          tags$tr(tags$td("PMIDs"), tags$td(ev$`Number of PMIDs`[1]), tags$td(ev$`Number of PMIDs`[2]))
        )
      ),
      div(style = "font-size: 12px; color: #888; margin-bottom: 8px;",
          sprintf("Clinical Association Score %.2f · TME %.3f · Literature %.3f · RRA P %s",
                  pr$clinical_score, pr$TME_score, pr$literature_score, signif(pr$RRA_score, 3)))
    )
  })
  
  output$download_evidence <- downloadHandler(
    filename = function() {
      pr <- selected_pair()
      paste0("TACTIC-GBM_evidence_", safe_file_part(pr$Drug1), "_", safe_file_part(pr$Drug2), ".csv")
    },
    content = function(file) {
      utils::write.csv(evidence_table(), file, row.names = FALSE, na = "")
    }
  )
  
  ## PMIDs by drug ------------------------------------------------------
  pmid_selected <- reactive({
    req(input$pmid_drugs)
    pmid_long %>% dplyr::filter(Drug %in% input$pmid_drugs)
  })
  
  output$pmid_preview <- DT::renderDT({
    validate(need(length(input$pmid_drugs) > 0, "Choose one or more drugs to preview their PubMed IDs."))
    show <- pmid_selected() %>%
      dplyr::mutate(PMID = sprintf('<a href="%s" target="_blank">%s</a>', `PubMed link`, PMID)) %>%
      dplyr::select(Drug, `Dominant cell group`, PMID)
    DT::datatable(
      show, rownames = FALSE, escape = FALSE,
      options = list(dom = "tip", pageLength = 5, scrollX = TRUE)
    )
  })
  
  output$download_pmids <- downloadHandler(
    filename = function() {
      d <- input$pmid_drugs
      tag <- if (length(d) == 1) safe_file_part(d) else paste0(length(d), "_drugs")
      paste0("TACTIC-GBM_PMIDs_", tag, ".csv")
    },
    content = function(file) {
      out <- if (length(input$pmid_drugs)) pmid_selected() else pmid_long[0, ]
      utils::write.csv(out, file, row.names = FALSE, na = "")
    }
  )
  
  ## Quick-pick buttons for the manuscript doublets
  lapply(seq_len(nrow(featured_doublets)), function(i) {
    observeEvent(input[[paste0("featured_doublet_", i)]], {
      d1 <- featured_doublets$app_drug1[i]
      d2 <- featured_doublets$app_drug2[i]
      if (!is.na(d1)) updateSelectInput(session, "drug1_pc1", selected = d1)
      if (!is.na(d2)) updateSelectInput(session, "drug2_pc1", selected = d2)
    }, ignoreInit = TRUE)
  })
  
  ## Tutorial navigation
  tutorial_links <- c(
    intro_to_tutorial  = "Tutorial",
    tut_go_top         = "Top Combinations",
    tut_go_doublet     = "Doublet Survival Analysis",
    tut_go_single_drug = "Single Drug Analysis",
    tut_go_heatmap     = "Single Cell Heatmap",
    tut_go_drug_gene   = "Drug-Gene",
    tut_go_top_dl      = "Top Combinations",
    intro_go_top         = "Top Combinations",
    intro_go_doublet     = "Doublet Survival Analysis",
    intro_go_single_drug = "Single Drug Analysis",
    intro_go_single_gene = "Single Gene Analysis",
    intro_go_heatmap     = "Single Cell Heatmap",
    intro_go_drug_gene   = "Drug-Gene",
    intro_go_contact     = "Contact",
    intro_btn_tutorial   = "Tutorial",
    intro_btn_top        = "Top Combinations"
  )
  lapply(names(tutorial_links), function(id) {
    observeEvent(input[[id]], nav_select("main_tabs", tutorial_links[[id]]), ignoreInit = TRUE)
  })
  
  # Single Gene Analysis links in the tutorial open the CTSB example
  show_ctsb <- function() {
    if ("CTSB" %in% unique_genes) {
      updateSelectizeInput(session, "gene", selected = "CTSB")
    } else {
      warning("CTSB is not in the Single Gene Analysis gene list.")
    }
    nav_select("main_tabs", "Single Gene Analysis")
  }
  observeEvent(input$tut_go_single_gene, show_ctsb(), ignoreInit = TRUE)
  
  ## Tutorial worked example: GW-559090 + CHEMBL1321401
  select_if_present <- function(input_id, wanted, choices) {
    hit <- match_drug_names(wanted, choices)
    if (!is.na(hit)) updateSelectInput(session, input_id, selected = hit)
    invisible(hit)
  }
  example_drug1 <- featured_doublets$drug1[1]
  example_drug2 <- featured_doublets$drug2[1]
  
  observeEvent(input$ex_doublet, {
    select_if_present("drug1_pc1", example_drug1, drug_seed_ls)
    select_if_present("drug2_pc1", example_drug2, unique(Doublet_all_with_coef$Drug2))
    nav_select("main_tabs", "Doublet Survival Analysis")
  })
  observeEvent(input$ex_single_drug, {
    select_if_present("drug", example_drug1, drug_list)
    nav_select("main_tabs", "Single Drug Analysis")
  })
  observeEvent(input$ex_heatmap, {
    select_if_present("Drug1", example_drug1, unique(merge_by_group$drug))
    nav_select("main_tabs", "Single Cell Heatmap")
  })
  observeEvent(input$ex_gene, {
    show_ctsb()
  })
  observeEvent(input$ex_drug_gene, {
    hit <- match_drug_names(example_drug1, drug_gene_tab_choices)
    if (!is.na(hit)) updateSelectizeInput(session, "drug_gene_selected", selected = hit)
    nav_select("main_tabs", "Drug-Gene")
  })
  
  
  output$pc1SurvivalPlot <- renderPlot({
    
    req(input$drug1_pc1, input$drug2_pc1)
    
    drug1 <- input$drug1_pc1
    drug2 <- input$drug2_pc1
    
    plots <- list()
    
    for (cohort_name in c("tcga", "cgga325")) {
      
      seq_final <- geneData()[[cohort_name]]
      Clinical_Data <- ClinicalData()[[cohort_name]]
      
      for (drug_name in c(drug1, drug2)) {
        
        p <- run_survival_plot_single_drug(
          seq_final = seq_final,
          Clinical_Data = Clinical_Data,
          drug_name = drug_name,
          cohort_name = toupper(cohort_name),
          cutoff_method = "Maxstat Cut-off"
        )
        
        p$plot <- p$plot + ggtitle(paste0(toupper(cohort_name), ": ", drug_name))
        
        plots[[paste(cohort_name, drug_name, sep = "_")]] <- p$plot
      }
    }
    
    # Row 1 = TCGA (drug1, drug2); Row 2 = CGGA325 (drug1, drug2)
    grid.arrange(
      plots[[paste("tcga",    drug1, sep = "_")]],
      plots[[paste("tcga",    drug2, sep = "_")]],
      plots[[paste("cgga325", drug1, sep = "_")]],
      plots[[paste("cgga325", drug2, sep = "_")]],
      ncol = 2, nrow = 2
    )
  })
  
  
  output$HeatmapPlot_by_type <- renderImage({
    req(input$Drug1)
    render_heatmap_image(
      plot_drug_heatmap_type,
      df = merge_by_type, drug_name = input$Drug1, gbm_datasets = gbm_datasets,
      cellwidth = 14, cellheight = 14,
      alt = "Cell-type heatmap"
    )
  }, deleteFile = TRUE)
  
  output$HeatmapPlot_by_group <- renderImage({
    req(input$Drug1)
    render_heatmap_image(
      plot_drug_heatmap_update,
      df = merge_by_group_5, drug_name = input$Drug1, gbm_datasets = gbm_datasets,
      cellwidth = 34, cellheight = 14,
      alt = "Cell-group heatmap"
    )
  }, deleteFile = TRUE)
  
  
  output$drugGeneTable <- DT::renderDT({
    
    gene_sets <- selected_drug_gene_sets()
    
    table_data <- tibble::tibble(
      Source = names(gene_sets),
      `Number of Target Genes` = lengths(gene_sets),
      `Target Genes` = vapply(
        gene_sets,
        function(genes) {
          if (length(genes) == 0) {
            "Not available"
          } else {
            paste(sort(genes), collapse = ", ")
          }
        },
        character(1)
      )
    )
    
    DT::datatable(
      table_data,
      rownames = FALSE,
      escape = TRUE,
      options = list(
        dom = "t",
        pageLength = 4,
        scrollX = TRUE,
        autoWidth = TRUE,
        columnDefs = list(
          list(width = "15%", targets = 0),
          list(width = "15%", targets = 1),
          list(width = "70%", targets = 2)
        )
      )
    )
  })
  
  output$drugGeneVenn1 <- renderPlot({
    print(plot_list[[1]])
  })
  
  
  output$drugGeneVenn2 <- renderPlot({
    print(plot_list[[2]])
  })
  
  
  output$drugGeneVenn3 <- renderPlot({
    print(plot_list[[3]])
  })
  
  
  output$drugGeneVenn4 <- renderPlot({
    print(plot_list[[4]])
  })
  
  output$GeneHeatmapPlot_by_group <- renderImage({
    req(input$Gene1)
    render_heatmap_image(
      plot_gene_heatmap_by_group,
      df = merge_by_group_gene, gene_name = input$Gene1, gbm_datasets = gbm_datasets,
      cellwidth = 34, cellheight = 14,
      alt = "Cell-group heatmap by gene"
    )
  }, deleteFile = TRUE)
  
  output$drug1_full_name_note <- renderUI({
    req(input$Drug1)
    if (nchar(input$Drug1) > 20) {
      div(
        style = "margin-top: 8px; text-align: center; color: #888; font-size: 12px; font-style: italic;",
        paste0("Full name: ", input$Drug1)
      )
    }
  })
  
  
  output$GeneHeatmapPlot_by_type <- renderImage({
    req(input$Gene1)
    render_heatmap_image(
      plot_gene_heatmap_by_type,
      df = merge_by_type_gene, gene_name = input$Gene1, gbm_datasets = gbm_datasets,
      cellwidth = 14, cellheight = 14,
      alt = "Cell-type heatmap by gene"
    )
  }, deleteFile = TRUE)
  
}

shinyApp(ui,server)
