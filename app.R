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


tcga_clinical_subtype = readRDS("data/tcga_clinical_subtype.rds")
cgga325_clinical_subtype = readRDS("data/CGGA_clinical_325_subtype.rds")
#cgga693_clinical_subtype = readRDS("data/CGGA_clinical_693_subtype.rds")

gbm_tcga_seq_final_geneanlysis <- readRDS("data/gbm_tcga_seq_final_geneanlysis_new.rds")
#CGGA_rna_seq_final_693_geneanlysis <- readRDS("data/CGGA_693_readcount_Dseq2.rds")
CGGA_rna_seq_final_325_geneanlysis <- readRDS("data/CGGA_325_readcount_Dseq2.rds")

####### Primary dataset for use ###########
gbm_tcga_seq_final_geneanlysis_primary = gbm_tcga_seq_final_geneanlysis %>% filter(patient_barcode %in% tcga_clinical_subtype$patient_barcode) 
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

plot_gene_heatmap_by_group <- function(df, gene_name, gbm_datasets, scale_rows = "none",
                                       cluster_rows = TRUE, cluster_cols = FALSE,
                                       na_col = "grey85", display_numbers = FALSE) {
  
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
    main = paste0(gene_name, " - % Expressing by Cell Group")
  )
}


plot_gene_heatmap_by_type <- function(df, gene_name, gbm_datasets, scale_rows = "none",
                                      cluster_rows = TRUE, cluster_cols = FALSE,
                                      na_col = "grey85", display_numbers = FALSE) {
  
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
    main = paste0(gene_name, " - % Expressing by Cell Type")
  )
}

plot_drug_heatmap_update <- function(df, drug_name, gbm_datasets, scale_rows = "none",
                                     cluster_rows = TRUE, cluster_cols = FALSE,
                                     na_col = "grey85", display_numbers = FALSE) {
  
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
    scale = scale_rows,                          # "none" -> raw % values, real color scale
    cluster_rows = FALSE,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    na_col = na_col,
    display_numbers = display_numbers,            # TRUE -> overlay actual % on each cell
    number_format = "%.1f",
    fontsize_col = 14,
    fontsize = 12,
    main = paste0(drug_name, " - % Expressing by Cell Group")
  )
}


plot_drug_heatmap_type <- function(df, drug_name, gbm_datasets, scale_rows = "none",
                                   cluster_rows = TRUE, cluster_cols = FALSE,
                                   na_col = "grey85", display_numbers = FALSE) {
  
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
    fontsize_row = 12,
    fontsize = 12,
    main = paste0(drug_name, " - % Expressing by Cell Type")
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
    "))
  ),
  
  div(class = "app-title", "TACTIC-GBM: Tumor microenvironment – Aware Combination Therapy Index and Characterization for glioblastoma"),
  div(class = "app-subtitle", "Interactive exploration of genes, single drugs, doublets, and single-cell patterns"),
  
  navset_card_tab(
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
          
          h5(strong("What you can do here")),
          
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
              strong("PC1 and Survival — "),
              "Pick any two drugs to jointly stratify patients by their combined PCA signatures",
              "(Low/High for each), and compare four-group survival outcomes across cohorts."
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
          
          p(
            style = "color: #888; font-size: 13px;",
            "Use the tabs above to navigate between analyses. Most tabs include selectable",
            "inputs (genes, drugs, cut-off methods) that update the plots and tables in real time."
          )
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
                options = list(placeholder = "Search for a drug...")
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
            selectInput("gene", "Gene", choices = unique_genes),
            plotOutput("survivalPlot", height = "500px")
          )
        )
      )
    ),
    
    nav_panel(
      "Single Drug Analysis",
      layout_column_wrap(
        width = 1,
        card(
          card_header("Single Drug Survival"),
          card_body(
            fluidRow(
              column(6, selectInput("drug", "Drug", choices = drug_list)),
              column(6, selectInput("cutoff_method", "Cut-off method",
                                    choices = c("Median Cut", "Maxstat Cut-off")))
            ),
            br(),
            fluidRow(
              column(6, card(card_header("TCGA Loadings"), card_body(DTOutput("pcLoading_tcga")))),
              column(6, card(card_header("CGGA325 Loadings"), card_body(DTOutput("pcLoading_cgga325"))))
            ),
            br(),
            card(
              card_header("Survival Curves"),
              card_body(plotOutput("drugSurvivalPlot", height = "520px"))
            )
          )
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
            p("Ranked candidate combinations with coefficients, p-values, and composite scores."),
            DTOutput("topCombTable")
          )
        )
      )
    ),
    
    nav_panel(
      "PC1 and Survival",
      layout_column_wrap(
        width = 1,
        card(
          card_header("Doublet Survival Analysis"),
          card_body(
            fluidRow(
              column(6, selectInput("drug1_pc1", "Drug 1", choices = drug_seed_ls)),
              column(6, selectInput("drug2_pc1", "Drug 2", choices = unique(Doublet_all_with_coef$Drug2)))
            ),
            br(),
            plotOutput("pc1SurvivalPlot", height = "520px")
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
            selectInput("Drug1", "Drug", choices = unique(merge_by_group$drug)),
            br(),
            fluidRow(
              column(6, card(card_header("Cell Group"), card_body(imageOutput("HeatmapPlot_by_group", height = "600px")))),
              column(6, card(card_header("Cell Type"), card_body(imageOutput("HeatmapPlot_by_type", height = "600px"))))
            )
          )
        ),
        
        card(
          card_header("Mechanistic Insights — By Gene"),
          card_body(
            selectInput("Gene1", "Gene", choices = unique(merge_by_group_gene$gene)),
            br(),
            fluidRow(
              column(6, card(card_header("Cell Group"), card_body(imageOutput("GeneHeatmapPlot_by_group", height = "600px")))),
              column(6, card(card_header("Cell Type"), card_body(imageOutput("GeneHeatmapPlot_by_type", height = "600px"))))
            )
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
  
  DrugData = reactive({
    selected_drug = input$drug
    
    drug_gene_ls = readRDS("data/drug_gene_list.rds")
    drug_list = readRDS("data/drug_list.rds")
    
    list(drug_gene_ls = drug_gene_ls,
         drug_list = drug_list,
         drug = selected_drug)
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
    
    plots <- list()
    drug_name <- input$drug
    cutoff_method <- input$cutoff_method
    
    for (cohort_name in c("tcga",  "cgga325")) {
      
      seq_final <- geneData()[[cohort_name]]
      Clinical_Data <- ClinicalData()[[cohort_name]]
      
      p <- run_survival_plot_single_drug(
        seq_final = seq_final,
        Clinical_Data = Clinical_Data,
        drug_name = drug_name,
        cohort_name = toupper(cohort_name),
        cutoff_method = cutoff_method
      )
      
      plots[[cohort_name]] <- p$plot
    }
    
    gridExtra::grid.arrange(
      grobs = plots,
      ncol = 2
    )
  })
  
  #drug_name = drug_list[1]
  
  output$pcLoading_tcga <- DT::renderDT({
    drug_name <- DrugData()$drug
    target_gene_ls <- DrugData()$drug_gene_ls[[drug_name]]
    
    pca_result_tcga <- data_PC(geneData()$tcga, target_gene_ls, drug_name)
    loading_tcga <- as.data.frame(pca_result_tcga$loading[[drug_name]]) %>%
      tibble::rownames_to_column("Gene") %>%
      head(10)
    
    datatable(
      loading_tcga,
      rownames = FALSE,
      options = list(dom = "t", pageLength = 10, scrollX = TRUE)
    ) %>%
      formatRound(columns = setdiff(names(loading_tcga), "Gene"), digits = 3)
  })
  

  
  output$pcLoading_cgga325 <- DT::renderDT({
    drug_name <- DrugData()$drug
    target_gene_ls <- DrugData()$drug_gene_ls[[drug_name]]
    
    pca_result_cgga325 <- data_PC(geneData()$cgga325, target_gene_ls, drug_name)
    loading_cgga325 <- as.data.frame(pca_result_cgga325$loading[[drug_name]]) %>%
      tibble::rownames_to_column("Gene") %>%
      head(10)
    
    datatable(
      loading_cgga325,
      rownames = FALSE,
      options = list(dom = "t", pageLength = 10, scrollX = TRUE)
    ) %>%
      formatRound(columns = setdiff(names(loading_cgga325), "Gene"), digits = 3)
  })
  
  output$topCombTable <- DT::renderDT({
    tab <- candidate_doublet %>%
      dplyr::arrange((RRA_score)) %>%
      dplyr::rename(
        `Drug 1` = Drug1,
        `Drug 1 Source` = drug_set_drug1,
        `Drug 2` = Drug2,
        `Drug 2 Source` = drug_set_drug2,
        `Clinical Association Score` = Fisher_combined_tcga_cgga325,
        `Tumor Microenvironment Score` = TME_score,
        `Literature Support Score` = literature_score,
        `Recipocal Rank Fusion Score` = RRF_score,
        `Robust Rank Aggregation Score` = RRA_score
      )
    
    datatable(
      tab,
      rownames = FALSE,
      filter = "top",
      extensions = c("Buttons", "Scroller"),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel"),
        autoWidth = TRUE
      )
    ) %>%
      formatRound(
        columns = names(tab)[sapply(tab, is.numeric)],
        digits = 3
      )
  })

  
  output$pc1SurvivalPlot = renderPlot({
    
    data_PC_doublet = function(data, target_gene_ls,drug_list){
      #Filter the gene expression data
      filter_data = data %>% 
        dplyr::select(any_of(unlist(target_gene_ls,use.names = FALSE)))
      
      loading = list()
      
      for (i in seq_along(drug_list)){
        pca = prcomp(filter_data[, colnames(filter_data) %in% target_gene_ls[[drug_list[i]]]], center = TRUE, scale = TRUE)
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
    
    drug_gene_ls = DrugData()$drug_gene_ls
    drug1 = input$drug1_pc1
    #drug1 = drug_seed_ls[1]
    drug2 = input$drug2_pc1
    #drug2 = drug_list[2]
    
    drug_list = c(drug1,drug2)
    
    plots <- list()
    for (cohort_name in c("tcga",  "cgga325")) {
      
      #corhot = "tcga"
      #seq_final <- gbm_tcga_seq_final_geneanlysis
      #Clinical_Data <- tcga_clinical
      seq_final <- geneData()[[cohort_name]]
      Clinical_Data <- ClinicalData()[[cohort_name]]
      #Clinical_Data <- tcga_clinical
      
      
      drug_target_gene_ls_doublet <- list()
      for (current_drug in drug_list) {
        # Assign the target genes to db_drug_target_gene_ls for this drug
        drug_target_gene_ls_doublet[[current_drug]] <- unique_drug_genes[[current_drug]]
      }
      
      filter_data <- add_column(data_PC_doublet(seq_final,drug_target_gene_ls_doublet,drug_list)$filter_data[,c("PC_1","PC_2")], seq_final["patient_barcode"]) 
      final_data <- left_join(filter_data, Clinical_Data, by="patient_barcode")
      
      cutpoint_1 <- median(final_data[["PC_1"]], na.rm = TRUE)
      cutpoint_2 <-median(final_data[["PC_2"]], na.rm = TRUE)
      
      df_sar_group <- final_data %>% 
        dplyr::select(os_time, os_status, PC_1, PC_2) %>% 
        dplyr::mutate(drug_1_group = if_else(PC_1 <= cutpoint_1, "Low", "High")) %>% 
        dplyr::mutate(drug_1_group = factor(drug_1_group, levels = c("Low", "High"))) %>% 
        
        dplyr::mutate(drug_2_group = if_else(PC_2 <= cutpoint_2, "Low", "High")) %>% 
        dplyr::mutate(drug_2_group = factor(drug_2_group, levels = c("Low", "High"))) 
      
      surv.model <- survfit(formula = Surv(os_time, os_status) ~ drug_1_group + drug_2_group,
                            data = df_sar_group)
      
      p = ggsurvplot(surv.model, risk.table = TRUE,
                     data = df_sar_group,
                     pval = TRUE,
                     font.tickslab = c(18),
                     ggtheme = theme_classic2(base_size=14), 
                     risk.table.fontsize = 4,  # Reduce font size for the risk table
                     risk.table.height = 0.4,  # Reduce the height of the risk table
                     tables.height = 0.3, xlab ="Survival time in Months",
                     legend.labs = c("Low,Low", "Low,High","High,Low","High,High"),
                     xlim = c(0, 60)
      )
      
      p$plot <- p$plot + ggtitle(toupper(cohort_name))
      
      plots[[cohort_name]] <- p$plot
      
    }
    
    grid.arrange(plots$tcga, plots$cgga325, ncol = 2)
    
  })
  
  
  output$HeatmapPlot_by_type <- renderImage({
    req(input$Drug1)
    
    outfile <- tempfile(fileext = ".png")
    
    png(outfile, width = 1400, height = 1100, res = 180)   # <- bumped up
    
    plot_drug_heatmap_type(
      df = merge_by_type,
      drug_name = input$Drug1,
      gbm_datasets = gbm_datasets
    )
    
    dev.off()
    
    list(
      src = outfile,
      contentType = "image/png",
      width = "100%",
      height = "600px",
      alt = "Cell-type heatmap"
    )
  }, deleteFile = TRUE)
  
  output$HeatmapPlot_by_group <- renderImage({
    
    req(input$Drug1)
    
    outfile <- tempfile(fileext = ".png")
    
    png(outfile, width = 1200, height = 900, res = 150)
    
    plot_drug_heatmap_update(
      df = merge_by_group_5,
      drug_name = input$Drug1,
      gbm_datasets = gbm_datasets
    )
    
    dev.off()
    
    list(
      src = outfile,
      contentType = "image/png",
      width = "100%",
      height = "600px",
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
    
    outfile <- tempfile(fileext = ".png")
    
    png(outfile, width = 1400, height = 1100, res = 180)
    
    plot_gene_heatmap_by_group(
      df = merge_by_group_gene,
      gene_name = input$Gene1,
      gbm_datasets = gbm_datasets
    )
    
    dev.off()
    
    list(
      src = outfile,
      contentType = "image/png",
      width = "100%",
      height = "600px",
      alt = "Cell-group heatmap by gene"
    )
  }, deleteFile = TRUE)
  
  
  output$GeneHeatmapPlot_by_type <- renderImage({
    req(input$Gene1)
    
    outfile <- tempfile(fileext = ".png")
    
    png(outfile, width = 1400, height = 1100, res = 180)
    
    plot_gene_heatmap_by_type(
      df = merge_by_type_gene,
      gene_name = input$Gene1,
      gbm_datasets = gbm_datasets
    )
    
    dev.off()
    
    list(
      src = outfile,
      contentType = "image/png",
      width = "100%",
      height = "600px",
      alt = "Cell-type heatmap by gene"
    )
  }, deleteFile = TRUE)
  
}

shinyApp(ui,server)


