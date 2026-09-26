# ============================================================
# 05_expression_preprocessing.R
# Clinical–Molecular Integration Project
# TCGA Pancreatic Adenocarcinoma (PAAD)
# ============================================================


# -----------------------------
# 1. Project paths
# -----------------------------

project_dir <- "E:/clinical-molecular-integration"

expression_file <- file.path(
  project_dir,
  "data",
  "data_mrna_seq_v2_rsem.txt"
)


# -----------------------------
# 2. Check input file
# -----------------------------

if (!file.exists(expression_file)) {
  stop("Expression file was not found.")
}

cat("Expression file found.\n\n")


# -----------------------------
# 3. Read expression data
# -----------------------------

expression_data <- read.delim(
  expression_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)


# -----------------------------
# 4. Identify annotation
#    and expression columns
# -----------------------------

gene_symbols <- expression_data$Hugo_Symbol

entrez_ids <- expression_data$Entrez_Gene_Id

expression_sample_ids <- colnames(
  expression_data
)[-c(1, 2)]


# -----------------------------
# 5. Basic preprocessing QC
# -----------------------------

cat("============================================\n")
cat("EXPRESSION PREPROCESSING QC\n")
cat("============================================\n\n")

cat(
  "Genes/features:",
  nrow(expression_data),
  "\n"
)

cat(
  "Samples:",
  length(expression_sample_ids),
  "\n"
)

cat(
  "Missing expression values:",
  sum(
    is.na(
      expression_data[, -(1:2)]
    )
  ),
  "\n"
)

cat(
  "Negative expression values:",
  sum(
    expression_data[, -(1:2)] < 0,
    na.rm = TRUE
  ),
  "\n"
)


# -----------------------------
# 6. Expression distribution
#    before transformation
# -----------------------------

expression_matrix <- as.matrix(
  expression_data[, -(1:2)]
)

expression_values <- as.numeric(
  expression_matrix
)


cat("\n============================================\n")
cat("EXPRESSION DISTRIBUTION\n")
cat("============================================\n\n")

cat(
  "Minimum:",
  min(expression_values, na.rm = TRUE),
  "\n"
)

cat(
  "1st Quartile:",
  quantile(
    expression_values,
    0.25,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Median:",
  median(
    expression_values,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Mean:",
  mean(
    expression_values,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "3rd Quartile:",
  quantile(
    expression_values,
    0.75,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Maximum:",
  max(expression_values, na.rm = TRUE),
  "\n"
)


# -----------------------------
# 7. Log2 transformation
# -----------------------------

expression_log2 <- log2(
  expression_matrix + 1
)


# -----------------------------
# 8. Post-transformation QC
# -----------------------------

expression_log2_values <- as.numeric(
  expression_log2
)

cat("\n============================================\n")
cat("POST-TRANSFORMATION QC\n")
cat("============================================\n\n")

cat(
  "Minimum:",
  min(expression_log2_values, na.rm = TRUE),
  "\n"
)

cat(
  "1st Quartile:",
  quantile(
    expression_log2_values,
    0.25,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Median:",
  median(
    expression_log2_values,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Mean:",
  mean(
    expression_log2_values,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "3rd Quartile:",
  quantile(
    expression_log2_values,
    0.75,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Maximum:",
  max(
    expression_log2_values,
    na.rm = TRUE
  ),
  "\n"
)


# -----------------------------
# 9. Gene-level expression QC
# -----------------------------

gene_medians <- apply(
  expression_log2,
  1,
  median,
  na.rm = TRUE
)


cat("\n============================================\n")
cat("GENE-LEVEL EXPRESSION QC\n")
cat("============================================\n\n")

cat(
  "Genes with median expression = 0:",
  sum(gene_medians == 0),
  "\n"
)

cat(
  "Genes with median expression < 1:",
  sum(gene_medians < 1),
  "\n"
)

cat(
  "Genes with median expression < 2:",
  sum(gene_medians < 2),
  "\n"
)

cat(
  "Genes with median expression < 3:",
  sum(gene_medians < 3),
  "\n"
)

cat(
  "Genes with median expression >= 5:",
  sum(gene_medians >= 5),
  "\n"
)


# -----------------------------
# 10. Gene expression prevalence
# -----------------------------

gene_detected_count <- rowSums(
  expression_matrix > 0,
  na.rm = TRUE
)

gene_detected_fraction <- (
  gene_detected_count /
    ncol(expression_matrix)
)


cat("\n============================================\n")
cat("GENE EXPRESSION PREVALENCE\n")
cat("============================================\n\n")

cat(
  "Genes detected in at least 10% of samples:",
  sum(gene_detected_fraction >= 0.10),
  "\n"
)

cat(
  "Genes detected in at least 25% of samples:",
  sum(gene_detected_fraction >= 0.25),
  "\n"
)

cat(
  "Genes detected in at least 50% of samples:",
  sum(gene_detected_fraction >= 0.50),
  "\n"
)

cat(
  "Genes detected in at least 75% of samples:",
  sum(gene_detected_fraction >= 0.75),
  "\n"
)

cat(
  "Genes detected in all samples:",
  sum(gene_detected_fraction == 1),
  "\n"
)


# -----------------------------
# 11. Gene filtering
# -----------------------------

minimum_samples <- ceiling(
  0.10 * ncol(expression_matrix)
)

gene_filter <- (
  gene_detected_count >= minimum_samples
)


cat("\n============================================\n")
cat("GENE FILTERING\n")
cat("============================================\n\n")

cat(
  "Total genes before filtering:",
  nrow(expression_matrix),
  "\n"
)

cat(
  "Minimum required samples:",
  minimum_samples,
  "\n"
)

cat(
  "Genes retained:",
  sum(gene_filter),
  "\n"
)

cat(
  "Genes removed:",
  sum(!gene_filter),
  "\n"
)


expression_log2_filtered <- expression_log2[
  gene_filter,
  ,
  drop = FALSE
]

gene_symbols_filtered <- gene_symbols[
  gene_filter
]

entrez_ids_filtered <- entrez_ids[
  gene_filter
]


cat(
  "Filtered expression matrix:",
  nrow(expression_log2_filtered),
  "genes x",
  ncol(expression_log2_filtered),
  "samples\n"
)


# -----------------------------
# 12. Post-filter gene ID QC
# -----------------------------

filtered_symbols_nonempty <- gene_symbols_filtered[
  !is.na(gene_symbols_filtered) &
    gene_symbols_filtered != ""
]

filtered_entrez_nonempty <- entrez_ids_filtered[
  !is.na(entrez_ids_filtered) &
    entrez_ids_filtered != ""
]


duplicate_filtered_symbols <- unique(
  filtered_symbols_nonempty[
    duplicated(filtered_symbols_nonempty)
  ]
)

duplicate_filtered_entrez <- unique(
  filtered_entrez_nonempty[
    duplicated(filtered_entrez_nonempty)
  ]
)


cat("\n============================================\n")
cat("POST-FILTER GENE ID QC\n")
cat("============================================\n\n")

cat(
  "Filtered genes:",
  nrow(expression_log2_filtered),
  "\n"
)

cat(
  "Missing/empty Hugo Symbols:",
  sum(
    is.na(gene_symbols_filtered) |
      gene_symbols_filtered == ""
  ),
  "\n"
)

cat(
  "Unique non-empty Hugo Symbols:",
  length(
    unique(filtered_symbols_nonempty)
  ),
  "\n"
)

cat(
  "Duplicated Hugo Symbols:",
  length(
    duplicate_filtered_symbols
  ),
  "\n"
)

cat(
  "Duplicated Entrez IDs:",
  length(
    duplicate_filtered_entrez
  ),
  "\n"
)

if (length(duplicate_filtered_symbols) > 0) {
  cat(
    "\nDuplicated Hugo Symbols:\n"
  )
  print(duplicate_filtered_symbols)
}


# -----------------------------
# 13. Final preprocessing QC
# -----------------------------

cat("\n============================================\n")
cat("FINAL PREPROCESSING QC\n")
cat("============================================\n\n")

cat(
  "Final genes:",
  nrow(expression_log2_filtered),
  "\n"
)

cat(
  "Final samples:",
  ncol(expression_log2_filtered),
  "\n"
)

cat(
  "Missing values:",
  sum(
    is.na(expression_log2_filtered)
  ),
  "\n"
)

cat(
  "Negative values:",
  sum(
    expression_log2_filtered < 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Sample IDs preserved:",
  identical(
    colnames(expression_log2_filtered),
    expression_sample_ids
  ),
  "\n"
)


# -----------------------------
# 14. Save preprocessing output
# -----------------------------

output_dir <- file.path(
  project_dir,
  "results",
  "preprocessed"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


preprocessed_expression <- list(
  expression = expression_log2_filtered,
  gene_symbols = gene_symbols_filtered,
  entrez_ids = entrez_ids_filtered,
  sample_ids = expression_sample_ids,
  filtering_threshold = minimum_samples,
  filtering_fraction = 0.10
)


saveRDS(
  preprocessed_expression,
  file.path(
    output_dir,
    "expression_preprocessed.rds"
  )
)


cat("\n============================================\n")
cat("PREPROCESSING OUTPUT\n")
cat("============================================\n\n")

cat(
  "Output file:",
  file.path(
    output_dir,
    "expression_preprocessed.rds"
  ),
  "\n"
)

cat(
  "File saved:",
  file.exists(
    file.path(
      output_dir,
      "expression_preprocessed.rds"
    )
  ),
  "\n"
)
