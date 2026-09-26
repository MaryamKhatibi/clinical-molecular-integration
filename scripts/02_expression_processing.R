# ============================================================
# 02_expression_processing.R
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
# 3. Read expression matrix
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
# 4. Basic dimensions
# -----------------------------

cat("============================================\n")
cat("EXPRESSION MATRIX\n")
cat("============================================\n\n")

cat("Number of rows:", nrow(expression_data), "\n")
cat("Number of columns:", ncol(expression_data), "\n\n")


# -----------------------------
# 5. Identify sample columns
# -----------------------------

expression_sample_ids <- colnames(expression_data)[-c(1, 2)]

cat("============================================\n")
cat("SAMPLE INFORMATION\n")
cat("============================================\n\n")

cat(
  "Number of expression samples:",
  length(expression_sample_ids),
  "\n\n"
)

cat("First 10 sample IDs:\n\n")

print(
  expression_sample_ids[
    1:min(10, length(expression_sample_ids))
  ]
)


# -----------------------------
# 6. Gene identifier QC
# -----------------------------

hugo_symbols <- expression_data$Hugo_Symbol
entrez_ids <- expression_data$Entrez_Gene_Id

valid_hugo <- !is.na(hugo_symbols) & hugo_symbols != ""

duplicate_symbols <- unique(
  hugo_symbols[
    duplicated(hugo_symbols) &
      !is.na(hugo_symbols) &
      hugo_symbols != ""
  ]
)
cat("\n============================================\n")
cat("GENE IDENTIFIER QC\n")
cat("============================================\n\n")

cat(
  "Missing Hugo Symbols:",
  sum(!valid_hugo),
  "\n"
)

cat(
  "Missing Entrez IDs:",
  sum(is.na(entrez_ids) | entrez_ids == ""),
  "\n"
)

cat(
  "Unique Hugo Symbols:",
  length(unique(hugo_symbols[valid_hugo])),
  "\n"
)

cat(
  "Duplicated Hugo Symbols:",
  length(duplicate_symbols),
  "\n"
)

if (length(duplicate_symbols) > 0) {
  cat("\nDuplicated symbols:\n")
  print(duplicate_symbols)
}


# -----------------------------
# 7. Missing expression values
# -----------------------------

expression_values <- expression_data[, -(1:2)]

cat("\n============================================\n")
cat("MISSING VALUES\n")
cat("============================================\n\n")

cat(
  "Total missing expression values:",
  sum(is.na(expression_values)),
  "\n"
)


# -----------------------------
# 8. Expression value summary
# -----------------------------

expression_matrix <- as.matrix(expression_values)

cat("\n============================================\n")
cat("EXPRESSION VALUES\n")
cat("============================================\n\n")

cat(
  "Minimum:",
  min(expression_matrix, na.rm = TRUE),
  "\n"
)

cat(
  "Maximum:",
  max(expression_matrix, na.rm = TRUE),
  "\n"
)

cat(
  "Median:",
  median(expression_matrix, na.rm = TRUE),
  "\n"
)


# -----------------------------
# 9. Save QC summary
# -----------------------------

qc_summary <- data.frame(
  metric = c(
    "Expression features",
    "Expression samples",
    "Missing Hugo Symbols",
    "Missing Entrez IDs",
    "Duplicated Hugo Symbols",
    "Missing expression values"
  ),
  value = c(
    nrow(expression_data),
    length(expression_sample_ids),
    sum(!valid_hugo),
    sum(is.na(entrez_ids) | entrez_ids == ""),
    length(duplicate_symbols),
    sum(is.na(expression_matrix))
  )
)

results_dir <- file.path(project_dir, "results", "tables")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

write.csv(
  qc_summary,
  file.path(results_dir, "expression_qc_summary.csv"),
  row.names = FALSE
)


# -----------------------------
# End
# -----------------------------

cat("\n============================================\n")
cat("Expression QC completed successfully.\n")
cat("QC summary saved.\n")
cat("============================================\n")