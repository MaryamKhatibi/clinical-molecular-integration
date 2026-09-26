# ============================================================
# 01_data_preparation.R
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

required_files <- c(
  clinical_patient_file,
  clinical_sample_file,
  expression_file
)

if (!all(file.exists(required_files))) {
  stop("One or more required input files were not found.")
}

cat("All required input files found.\n\n")


# -----------------------------
# 3. Function for cBioPortal
#    clinical files
# -----------------------------

read_cbio_clinical <- function(file) {
  
  read.delim(
    file,
    sep = "\t",
    header = TRUE,
    skip = 4,
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}


# -----------------------------
# 4. Read clinical data
# -----------------------------

clinical_patient <- read_cbio_clinical(
  clinical_patient_file
)

clinical_sample <- read_cbio_clinical(
  clinical_sample_file
)


# -----------------------------
# 5. Read expression header only
# -----------------------------
# We do NOT load the entire expression
# matrix yet. We only need the sample IDs
# for the initial cohort check.

expression_header <- strsplit(
  readLines(expression_file, n = 1),
  "\t",
  fixed = TRUE
)[[1]]

expression_sample_ids <- expression_header[-c(1, 2)]


# -----------------------------
# 6. Basic dimensions
# -----------------------------

n_patients <- length(
  unique(clinical_patient$PATIENT_ID)
)

n_clinical_samples <- length(
  unique(clinical_sample$SAMPLE_ID)
)

n_expression_samples <- length(
  unique(expression_sample_ids)
)


# -----------------------------
# 7. Match clinical samples
#    with expression samples
# -----------------------------

common_samples <- intersect(
  clinical_sample$SAMPLE_ID,
  expression_sample_ids
)

n_common_samples <- length(
  unique(common_samples)
)


# -----------------------------
# 8. Identify integrated
#    patients
# -----------------------------

integrated_sample_data <- clinical_sample[
  clinical_sample$SAMPLE_ID %in% common_samples,
]

integrated_patient_ids <- unique(
  integrated_sample_data$PATIENT_ID
)

integrated_patient_data <- clinical_patient[
  clinical_patient$PATIENT_ID %in% integrated_patient_ids,
]

n_integrated_patients <- length(
  unique(integrated_patient_data$PATIENT_ID)
)


# -----------------------------
# 9. Overall Survival data
# -----------------------------

os_status_table <- table(
  integrated_patient_data$OS_STATUS,
  useNA = "ifany"
)

os_complete <- (
  !is.na(integrated_patient_data$OS_MONTHS) &
    is.finite(integrated_patient_data$OS_MONTHS) &
    integrated_patient_data$OS_MONTHS >= 0 &
    !is.na(integrated_patient_data$OS_STATUS) &
    nzchar(integrated_patient_data$OS_STATUS)
)

n_complete_os <- sum(os_complete)


# -----------------------------
# 10. Print cohort summary
# -----------------------------

cat("============================================\n")
cat("COHORT SUMMARY\n")
cat("============================================\n\n")

cat(
  "Unique patients in clinical data:",
  n_patients,
  "\n"
)

cat(
  "Unique samples in clinical sample data:",
  n_clinical_samples,
  "\n"
)

cat(
  "Samples in expression matrix:",
  n_expression_samples,
  "\n"
)

cat(
  "Samples shared between clinical and expression data:",
  n_common_samples,
  "\n"
)

cat(
  "Patients after clinical-expression integration:",
  n_integrated_patients,
  "\n"
)

cat(
  "Patients with complete OS information:",
  n_complete_os,
  "\n\n"
)


# -----------------------------
# 11. OS status distribution
# -----------------------------

cat("============================================\n")
cat("OVERALL SURVIVAL STATUS\n")
cat("============================================\n\n")

print(os_status_table)


# -----------------------------
# 12. Save summary tables
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

cohort_summary <- data.frame(
  metric = c(
    "clinical_patients",
    "clinical_samples",
    "expression_samples",
    "common_samples",
    "integrated_patients",
    "complete_overall_survival"
  ),
  value = c(
    n_patients,
    n_clinical_samples,
    n_expression_samples,
    n_common_samples,
    n_integrated_patients,
    n_complete_os
  )
)

write.csv(
  cohort_summary,
  file.path(
    results_dir,
    "cohort_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(os_status_table),
  file.path(
    results_dir,
    "os_status_summary.csv"
  ),
  row.names = FALSE
)


cat("\nSummary tables saved to:\n")
cat(results_dir, "\n")

# -----------------------------
# Identify samples without
# matching expression data
# -----------------------------

missing_expression_samples <- setdiff(
  clinical_sample$SAMPLE_ID,
  expression_sample_ids
)

missing_expression_data <- clinical_sample[
  clinical_sample$SAMPLE_ID %in% missing_expression_samples,
  c(
    "PATIENT_ID",
    "SAMPLE_ID",
    "SAMPLE_TYPE",
    "CANCER_TYPE",
    "CANCER_TYPE_DETAILED",
    "TUMOR_TISSUE_SITE"
  )
]

print(missing_expression_data)

# -----------------------------
# Check whether the 7 patients
# have any other expression
# samples
# -----------------------------

missing_patient_ids <- unique(
  missing_expression_data$PATIENT_ID
)

alternative_expression_samples <- clinical_sample[
  clinical_sample$PATIENT_ID %in% missing_patient_ids &
    clinical_sample$SAMPLE_ID %in% expression_sample_ids,
  c(
    "PATIENT_ID",
    "SAMPLE_ID",
    "SAMPLE_TYPE"
  )
]

print(alternative_expression_samples)