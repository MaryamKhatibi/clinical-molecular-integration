# ============================================================
# 03_build_analysis_dataset.R
# Clinical–Molecular Integration Project
# TCGA Pancreatic Adenocarcinoma (PAAD)
# ============================================================


# -----------------------------
# 1. Project paths
# -----------------------------

project_dir <- "E:/clinical-molecular-integration"

clinical_patient_file <- file.path(
  project_dir,
  "data",
  "data_clinical_patient.txt"
)

clinical_sample_file <- file.path(
  project_dir,
  "data",
  "data_clinical_sample.txt"
)

expression_file <- file.path(
  project_dir,
  "data",
  "data_mrna_seq_v2_rsem.txt"
)


# -----------------------------
# 2. Check input files
# -----------------------------

if (!file.exists(clinical_patient_file)) {
  stop("Clinical patient file was not found.")
}

if (!file.exists(clinical_sample_file)) {
  stop("Clinical sample file was not found.")
}

if (!file.exists(expression_file)) {
  stop("Expression file was not found.")
}

cat("Input files found.\n\n")


# -----------------------------
# 3. Read clinical data
# -----------------------------

clinical_patient <- read.delim(
  clinical_patient_file,
  sep = "\t",
  header = TRUE,
  skip = 4,
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

clinical_sample <- read.delim(
  clinical_sample_file,
  sep = "\t",
  header = TRUE,
  skip = 4,
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)


# -----------------------------
# 4. Select relevant variables
# -----------------------------

patient_variables <- clinical_patient[
  ,
  c(
    "PATIENT_ID",
    "AGE",
    "SEX",
    "AJCC_PATHOLOGIC_TUMOR_STAGE",
    "PATH_T_STAGE",
    "PATH_N_STAGE",
    "PATH_M_STAGE",
    "OS_STATUS",
    "OS_MONTHS"
  )
]

sample_variables <- clinical_sample[
  ,
  c(
    "PATIENT_ID",
    "SAMPLE_ID",
    "SAMPLE_TYPE",
    "TUMOR_TISSUE_SITE",
    "TUMOR_TYPE"
  )
]


# -----------------------------
# 5. Merge patient and sample data
# -----------------------------

clinical_data <- merge(
  sample_variables,
  patient_variables,
  by = "PATIENT_ID"
)


# -----------------------------
# 6. Identify expression samples
# -----------------------------

expression_header <- read.delim(
  expression_file,
  sep = "\t",
  header = TRUE,
  nrows = 1,
  check.names = FALSE
)

expression_sample_ids <- colnames(expression_header)[-c(1, 2)]


# -----------------------------
# 7. Keep samples with
#    expression data
# -----------------------------

analysis_clinical <- clinical_data[
  clinical_data$SAMPLE_ID %in% expression_sample_ids,
]


# -----------------------------
# 8. Integration checks
# -----------------------------

cat("============================================\n")
cat("CLINICAL–EXPRESSION INTEGRATION\n")
cat("============================================\n\n")

cat(
  "Clinical samples:",
  nrow(clinical_data),
  "\n"
)

cat(
  "Samples with expression data:",
  nrow(analysis_clinical),
  "\n"
)

cat(
  "Unique patients in analysis:",
  length(unique(analysis_clinical$PATIENT_ID)),
  "\n"
)

cat(
  "Complete OS information:",
  sum(
    !is.na(analysis_clinical$OS_MONTHS) &
      !is.na(analysis_clinical$OS_STATUS)
  ),
  "\n"
)


# -----------------------------
# 9. Select final variables
# -----------------------------

analysis_clinical <- analysis_clinical[
  ,
  c(
    "PATIENT_ID",
    "SAMPLE_ID",
    "SAMPLE_TYPE",
    "TUMOR_TISSUE_SITE",
    "TUMOR_TYPE",
    "AGE",
    "SEX",
    "AJCC_PATHOLOGIC_TUMOR_STAGE",
    "PATH_T_STAGE",
    "PATH_N_STAGE",
    "PATH_M_STAGE",
    "OS_STATUS",
    "OS_MONTHS"
  )
]


# -----------------------------
# 10. Clinical variable QC
# -----------------------------

clinical_variables <- c(
  "AGE",
  "SEX",
  "AJCC_PATHOLOGIC_TUMOR_STAGE",
  "PATH_T_STAGE",
  "PATH_N_STAGE",
  "PATH_M_STAGE",
  "OS_STATUS",
  "OS_MONTHS"
)

cat("\n============================================\n")
cat("CLINICAL VARIABLE QC\n")
cat("============================================\n\n")

cat("Missing values:\n\n")

for (variable in clinical_variables) {
  
  cat(
    variable,
    ":",
    sum(
      is.na(analysis_clinical[[variable]]) |
        analysis_clinical[[variable]] == ""
    ),
    "\n"
  )
}


# -----------------------------
# 11. Save final clinical dataset
# -----------------------------

results_dir <- file.path(
  project_dir,
  "results",
  "tables"
)

if (!dir.exists(results_dir)) {
  dir.create(
    results_dir,
    recursive = TRUE
  )
}

write.csv(
  analysis_clinical,
  file.path(
    results_dir,
    "analysis_clinical.csv"
  ),
  row.names = FALSE
)


# -----------------------------
# 12. Final check
# -----------------------------

cat("\n============================================\n")
cat("FINAL CLINICAL DATASET\n")
cat("============================================\n\n")

cat(
  "Rows:",
  nrow(analysis_clinical),
  "\n"
)

cat(
  "Columns:",
  ncol(analysis_clinical),
  "\n"
)

cat("\nClinical dataset saved successfully.\n")