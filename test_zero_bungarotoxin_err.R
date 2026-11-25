library(dplyr)

# Load the same objects app.R loads
CGGA_rna_seq_final_325_geneanlysis <- readRDS("data/CGGA_325_readcount_Dseq2.rds")
drug_gene_ls <- readRDS("data/drug_gene_list.rds")
drug_list    <- readRDS("data/drug_list.rds")

# Same drug/gene panel that triggers the crash
drug_name <- drug_list[grepl("Bungarotoxin", drug_list, ignore.case = TRUE)][1]
target_gene_ls <- drug_gene_ls[[drug_name]]

cat("Drug:", drug_name, "\n")
cat("Target genes:", paste(unlist(target_gene_ls), collapse = ", "), "\n\n")

# Same column selection data_PC() does
filter_data <- CGGA_rna_seq_final_325_geneanlysis %>%
  select(any_of(unlist(target_gene_ls, use.names = FALSE)))

cat("Columns actually present in this cohort:", paste(colnames(filter_data), collapse = ", "), "\n\n")

# The number that matters - print at full precision
col_sd <- sapply(filter_data, sd, na.rm = TRUE)
print(col_sd, digits = 20)


