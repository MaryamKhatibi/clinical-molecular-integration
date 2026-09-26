# ============================================================
# 04_expression_clinical_integration.R
# Clinical–Molecular Integration Project
# TCGA Pancreatic Adenocarcinoma (PAAD)
# ============================================================


# -----------------------------
# 1. Project paths
# -----------------------------

project_dir <- "E:/clinical-molecular-integration"

clinical_file <- file.path(
  project_dir,
  "results",
  "tables",
  "analysis_clinical.csv"
)

expression_file <- file.path(
  project_dir,
  "data",
  "data_mrna_seq_v2_rsem.txt"
)


# -----------------------------
# 2. Check input files
# -----------------------------

if (!file.exists(clinical_file)) {
  stop("Clinical analysis file was not found.")
}

if (!file.exists(expression_file)) {
  stop("Expression file was not found.")
}

cat("Input files found.\n\n")


# -----------------------------
# 3. Read clinical dataset
# -----------------------------

analysis_clinical <- read.csv(
  clinical_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Clinical dataset loaded.\n")
cat(
  "Rows:",
  nrow(analysis_clinical),
  "\n"
)

cat(
  "Columns:",
  ncol(analysis_clinical),
  "\n\n"
)


# -----------------------------
# 4. Read expression matrix
# -----------------------------

expression_data <- read.delim(
  expression_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

cat("Expression dataset loaded.\n")
cat(
  "Rows:",
  nrow(expression_data),
  "\n"
)

cat(
  "Columns:",
  ncol(expression_data),
  "\n\n"
)

# -----------------------------
# 5. Identify sample IDs
# -----------------------------

clinical_sample_ids <- analysis_clinical$SAMPLE_ID

expression_sample_ids <- colnames(expression_data)[-c(1, 2)]


# -----------------------------
# 6. Sample ID QC
# -----------------------------

cat("\n============================================\n")
cat("SAMPLE ID QC\n")
cat("============================================\n\n")

cat(
  "Clinical sample IDs:",
  length(clinical_sample_ids),
  "\n"
)

cat(
  "Expression sample IDs:",
  length(expression_sample_ids),
  "\n\n"
)


# Check duplicate IDs

cat(
  "Duplicated clinical sample IDs:",
  sum(duplicated(clinical_sample_ids)),
  "\n"
)

cat(
  "Duplicated expression sample IDs:",
  sum(duplicated(expression_sample_ids)),
  "\n\n"
)


# Find samples present in clinical data
# but missing from expression data

clinical_not_expression <- setdiff(
  clinical_sample_ids,
  expression_sample_ids
)


# Find samples present in expression data
# but missing from clinical data

expression_not_clinical <- setdiff(
  expression_sample_ids,
  clinical_sample_ids
)


cat(
  "Clinical samples not found in expression:",
  length(clinical_not_expression),
  "\n"
)

cat(
  "Expression samples not found in clinical:",
  length(expression_not_clinical),
  "\n\n"
)


# -----------------------------
# 7. Exact sample matching
# -----------------------------

if (
  length(clinical_not_expression) == 0 &&
  length(expression_not_clinical) == 0
) {
  
  cat("All sample IDs match exactly.\n")
  
} else {
  
  cat("Sample ID mismatch detected.\n")
  
}

# -----------------------------
# 8. Align expression samples
#    to clinical sample order
# -----------------------------

expression_sample_ids <- colnames(expression_data)[-c(1, 2)]

sample_order <- match(
  clinical_sample_ids,
  expression_sample_ids
)

if (any(is.na(sample_order))) {
  stop("Some clinical samples could not be matched to expression data.")
}


# Add 2 because expression_sample_ids
# starts after the first two annotation columns

expression_column_order <- sample_order + 2

expression_data_aligned <- expression_data[
  ,
  c(1, 2, expression_column_order)
]

aligned_sample_ids <- colnames(
  expression_data_aligned
)[-c(1, 2)]


# -----------------------------
# 9. Verify sample order
# -----------------------------

cat("\n============================================\n")
cat("SAMPLE ORDER QC\n")
cat("============================================\n\n")

cat(
  "Clinical samples:",
  length(clinical_sample_ids),
  "\n"
)

cat(
  "Aligned expression samples:",
  length(aligned_sample_ids),
  "\n"
)

cat(
  "Samples in identical order:",
  identical(
    clinical_sample_ids,
    aligned_sample_ids
  ),
  "\n"
)

if (
  !identical(
    clinical_sample_ids,
    aligned_sample_ids
  )
) {
  stop("Sample order does not match.")
}

cat("\nSample order aligned successfully.\n")

# -----------------------------
# 10. Final integration QC
# -----------------------------

cat("\n============================================\n")
cat("FINAL INTEGRATION QC\n")
cat("============================================\n\n")

cat(
  "Clinical rows:",
  nrow(analysis_clinical),
  "\n"
)

cat(
  "Expression genes/features:",
  nrow(expression_data_aligned),
  "\n"
)

cat(
  "Expression samples:",
  ncol(expression_data_aligned) - 2,
  "\n"
)

cat(
  "Clinical and expression sample counts match:",
  nrow(analysis_clinical) ==
    (ncol(expression_data_aligned) - 2),
  "\n"
)

cat(
  "Sample IDs and order match:",
  identical(
    analysis_clinical$SAMPLE_ID,
    colnames(expression_data_aligned)[-c(1, 2)]
  ),
  "\n"
)