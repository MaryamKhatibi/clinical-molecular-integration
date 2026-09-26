# ============================================================
# 06_expression_clinical_association.R
# Clinical–Molecular Integration Project
# TCGA Pancreatic Adenocarcinoma (PAAD)
# ============================================================

# ============================================================
# 1. PROJECT SETUP
# ============================================================

project_dir <- "E:/clinical-molecular-integration"

expression_file <- file.path(
  project_dir, "results", "preprocessed",
  "expression_preprocessed.rds"
)

clinical_file <- file.path(
  project_dir, "results", "tables",
  "analysis_clinical.csv"
)

results_dir <- file.path(
  project_dir, "results", "tables"
)

figures_dir <- file.path(
  project_dir, "figures"
)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is not installed.")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is not installed.")
}

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("Package 'ggrepel' is not installed.")
}

library(survival)
library(ggplot2)
library(ggrepel)

if (!file.exists(expression_file)) {
  stop("Preprocessed expression file was not found.")
}

if (!file.exists(clinical_file)) {
  stop("Clinical analysis file was not found.")
}

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}


# ============================================================
# 2. LOAD AND VALIDATE INPUT DATA
# ============================================================

preprocessed_expression <- readRDS(expression_file)

expression_matrix <- preprocessed_expression$expression
gene_symbols <- preprocessed_expression$gene_symbols
entrez_ids <- preprocessed_expression$entrez_ids
expression_sample_ids <- preprocessed_expression$sample_ids

analysis_clinical <- read.csv(
  clinical_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

clinical_sample_ids <- analysis_clinical$SAMPLE_ID

if (length(expression_sample_ids) != length(clinical_sample_ids)) {
  stop("Expression and clinical sample counts do not match.")
}

if (!identical(expression_sample_ids, clinical_sample_ids)) {
  stop("Clinical and expression sample IDs are not in identical order.")
}

n_genes <- nrow(expression_matrix)
n_samples <- ncol(expression_matrix)

cat("Input data loaded successfully.\n")
cat("Genes:", n_genes, "\n")
cat("Samples:", n_samples, "\n")


# ============================================================
# 3. CLINICAL VARIABLE PREPARATION
# ============================================================

analysis_clinical$AGE <- as.numeric(analysis_clinical$AGE)
analysis_clinical$SEX <- factor(analysis_clinical$SEX)

if (any(is.na(analysis_clinical$AGE))) {
  stop("AGE contains missing values.")
}

if (any(is.na(analysis_clinical$SEX)) ||
    any(analysis_clinical$SEX == "")) {
  stop("SEX contains missing or empty values.")
}

if (nlevels(analysis_clinical$SEX) != 2) {
  stop("SEX must contain exactly two groups.")
}

cat("\nClinical variables:\n")
cat("AGE range:",
    min(analysis_clinical$AGE),
    "to",
    max(analysis_clinical$AGE),
    "\n")
cat("SEX distribution:\n")
print(table(analysis_clinical$SEX))


# ============================================================
# 4. GENE EXPRESSION vs AGE
#    Spearman correlation + BH correction
# ============================================================

age <- analysis_clinical$AGE

rho_values <- numeric(n_genes)
p_values_age <- numeric(n_genes)

for (i in seq_len(n_genes)) {

  test_result <- cor.test(
    expression_matrix[i, ],
    age,
    method = "spearman",
    exact = FALSE
  )

  rho_values[i] <- unname(test_result$estimate)
  p_values_age[i] <- test_result$p.value
}

age_results <- data.frame(
  gene_symbol = gene_symbols,
  entrez_id = entrez_ids,
  rho = rho_values,
  p_value = p_values_age,
  FDR = p.adjust(p_values_age, method = "BH"),
  stringsAsFactors = FALSE
)

significant_age <- age_results[
  age_results$FDR < 0.05,
]

significant_age <- significant_age[
  order(significant_age$FDR),
]

write.csv(
  age_results,
  file.path(results_dir, "age_expression_association.csv"),
  row.names = FALSE
)

write.csv(
  significant_age,
  file.path(results_dir, "age_expression_association_FDR05.csv"),
  row.names = FALSE
)

cat("\nAge association:\n")
cat("Genes tested:", nrow(age_results), "\n")
cat("Nominal p < 0.05:",
    sum(age_results$p_value < 0.05, na.rm = TRUE),
    "\n")
cat("FDR < 0.10:",
    sum(age_results$FDR < 0.10, na.rm = TRUE),
    "\n")
cat("FDR < 0.05:",
    sum(age_results$FDR < 0.05, na.rm = TRUE),
    "\n")


# ============================================================
# 5. GENE EXPRESSION vs SEX
#    Wilcoxon rank-sum + rank-biserial effect size + BH
# ============================================================

sex <- analysis_clinical$SEX
group_1 <- levels(sex)[1]
group_2 <- levels(sex)[2]

p_values_sex <- numeric(n_genes)
effect_size_sex <- numeric(n_genes)

for (i in seq_len(n_genes)) {

  expression_values <- expression_matrix[i, ]

  group_1_values <- expression_values[sex == group_1]
  group_2_values <- expression_values[sex == group_2]

  test_result <- wilcox.test(
    group_1_values,
    group_2_values,
    exact = FALSE
  )

  p_values_sex[i] <- test_result$p.value

  n1 <- length(group_1_values)
  n2 <- length(group_2_values)

  combined_values <- c(
    group_1_values,
    group_2_values
  )

  ranks <- rank(
    combined_values,
    ties.method = "average"
  )

  rank_sum_1 <- sum(ranks[seq_len(n1)])
  U1 <- rank_sum_1 - (n1 * (n1 + 1) / 2)

  effect_size_sex[i] <- (
    2 * U1 / (n1 * n2)
  ) - 1
}

sex_results <- data.frame(
  gene_symbol = gene_symbols,
  entrez_id = entrez_ids,
  effect_size = effect_size_sex,
  p_value = p_values_sex,
  FDR = p.adjust(p_values_sex, method = "BH"),
  stringsAsFactors = FALSE
)

significant_sex <- sex_results[
  sex_results$FDR < 0.05,
]

significant_sex <- significant_sex[
  order(significant_sex$FDR),
]

write.csv(
  sex_results,
  file.path(results_dir, "sex_expression_association.csv"),
  row.names = FALSE
)

write.csv(
  significant_sex,
  file.path(results_dir, "sex_expression_association_FDR05.csv"),
  row.names = FALSE
)

cat("\nSex association:\n")
cat("Genes tested:", nrow(sex_results), "\n")
cat("Nominal p < 0.05:",
    sum(sex_results$p_value < 0.05, na.rm = TRUE),
    "\n")
cat("FDR < 0.10:",
    sum(sex_results$FDR < 0.10, na.rm = TRUE),
    "\n")
cat("FDR < 0.05:",
    sum(sex_results$FDR < 0.05, na.rm = TRUE),
    "\n")


# ============================================================
# 6. OVERALL SURVIVAL
#    Event coding + Kaplan–Meier estimate
# ============================================================

os_time <- analysis_clinical$OS_MONTHS
os_status <- analysis_clinical$OS_STATUS

event <- ifelse(
  os_status == "1:DECEASED",
  1,
  ifelse(os_status == "0:LIVING", 0, NA)
)

if (any(is.na(os_time)) ||
    any(is.na(event)) ||
    any(os_time < 0, na.rm = TRUE)) {
  stop("Survival data failed basic QC.")
}

cat("\nOverall survival:\n")
cat("N:", n_samples, "\n")
cat("Events:", sum(event == 1), "\n")
cat("Censored:", sum(event == 0), "\n")
cat("Median OS:",
    median(os_time),
    "months\n")

km_fit <- survfit(
  Surv(os_time, event) ~ 1
)

km_summary <- summary(km_fit)

median_os <- km_summary$table["median"]

cat("KM median OS:", median_os, "months\n")


# ============================================================
# 7. UNADJUSTED GENE-WISE COX REGRESSION
#    Exploratory reference analysis
# ============================================================

cox_p_values <- numeric(n_genes)
cox_hr <- numeric(n_genes)
cox_ci_lower <- numeric(n_genes)
cox_ci_upper <- numeric(n_genes)

for (i in seq_len(n_genes)) {

  cox_data <- data.frame(
    OS_MONTHS = os_time,
    event = event,
    expression = as.numeric(expression_matrix[i, ])
  )

  cox_model <- tryCatch(
    coxph(
      Surv(OS_MONTHS, event) ~ expression,
      data = cox_data
    ),
    error = function(e) NULL
  )

  if (is.null(cox_model)) {
    cox_p_values[i] <- NA
    cox_hr[i] <- NA
    cox_ci_lower[i] <- NA
    cox_ci_upper[i] <- NA
    next
  }

  model_summary <- summary(cox_model)

  cox_p_values[i] <- model_summary$coefficients[
    "expression", "Pr(>|z|)"
  ]

  cox_hr[i] <- model_summary$coefficients[
    "expression", "exp(coef)"
  ]

  cox_ci_lower[i] <- model_summary$conf.int[
    "expression", "lower .95"
  ]

  cox_ci_upper[i] <- model_summary$conf.int[
    "expression", "upper .95"
  ]
}

cox_results <- data.frame(
  gene_symbol = gene_symbols,
  entrez_id = entrez_ids,
  HR = cox_hr,
  CI_lower = cox_ci_lower,
  CI_upper = cox_ci_upper,
  p_value = cox_p_values,
  FDR = p.adjust(cox_p_values, method = "BH"),
  stringsAsFactors = FALSE
)

cox_results_sorted <- cox_results[
  order(cox_results$FDR),
]

significant_cox <- cox_results_sorted[
  cox_results_sorted$FDR < 0.05,
]

write.csv(
  cox_results,
  file.path(results_dir, "cox_survival_association.csv"),
  row.names = FALSE
)

write.csv(
  significant_cox,
  file.path(results_dir, "cox_survival_association_FDR05.csv"),
  row.names = FALSE
)

cat("\nUnadjusted Cox:\n")
cat("Genes tested:", nrow(cox_results), "\n")
cat("Successful models:",
    sum(!is.na(cox_results$p_value)),
    "\n")
cat("Failed models:",
    sum(is.na(cox_results$p_value)),
    "\n")
cat("FDR < 0.05:",
    sum(cox_results$FDR < 0.05, na.rm = TRUE),
    "\n")


# ============================================================
# 8. ADJUSTED COX REGRESSION
#    Expression + age + sex
# ============================================================

analysis_clinical$OS_event <- as.numeric(event)

adjusted_cox_results <- vector(
  "list",
  n_genes
)

for (i in seq_len(n_genes)) {

  cox_data <- data.frame(
    OS_MONTHS = analysis_clinical$OS_MONTHS,
    OS_event = analysis_clinical$OS_event,
    AGE = analysis_clinical$AGE,
    SEX = analysis_clinical$SEX,
    expression = as.numeric(expression_matrix[i, ])
  )

  fit <- tryCatch(
    coxph(
      Surv(OS_MONTHS, OS_event) ~
        expression + AGE + SEX,
      data = cox_data
    ),
    error = function(e) NULL
  )

  if (!is.null(fit)) {

    model_summary <- summary(fit)

    expression_row <- which(
      rownames(model_summary$coefficients) == "expression"
    )

    adjusted_cox_results[[i]] <- data.frame(
      gene_symbol = gene_symbols[i],
      entrez_id = entrez_ids[i],
      HR = model_summary$coefficients[
        expression_row, "exp(coef)"
      ],
      CI_lower = model_summary$conf.int[
        expression_row, "lower .95"
      ],
      CI_upper = model_summary$conf.int[
        expression_row, "upper .95"
      ],
      p_value = model_summary$coefficients[
        expression_row, "Pr(>|z|)"
      ]
    )
  }
}

adjusted_cox_results <- do.call(
  rbind,
  adjusted_cox_results
)

adjusted_cox_results$FDR <- p.adjust(
  adjusted_cox_results$p_value,
  method = "BH"
)

adjusted_cox_results <- adjusted_cox_results[
  order(adjusted_cox_results$FDR),
]

write.csv(
  adjusted_cox_results,
  file.path(
    results_dir,
    "cox_survival_association_adjusted_age_sex.csv"
  ),
  row.names = FALSE
)

cat("\nAdjusted Cox:\n")
cat("Genes successfully analyzed:",
    nrow(adjusted_cox_results),
    "\n")
cat("Nominal p < 0.05:",
    sum(adjusted_cox_results$p_value < 0.05),
    "\n")
cat("FDR < 0.10:",
    sum(adjusted_cox_results$FDR < 0.10),
    "\n")
cat("FDR < 0.05:",
    sum(adjusted_cox_results$FDR < 0.05),
    "\n")


# ============================================================
# 9. ADJUSTED COX CANDIDATE GENES
#    FDR < 0.05 and |log2(HR)| >= 0.5
# ============================================================

candidate_adjusted_cox <- subset(
  adjusted_cox_results,
  FDR < 0.05 &
    abs(log(HR, base = 2)) >= 0.5
)

candidate_adjusted_cox <- candidate_adjusted_cox[
  order(candidate_adjusted_cox$FDR),
]

write.csv(
  candidate_adjusted_cox,
  file.path(
    results_dir,
    "candidate_survival_genes_adjusted_age_sex.csv"
  ),
  row.names = FALSE
)

cat("\nAdjusted Cox candidates:\n")
cat("Candidate genes:",
    nrow(candidate_adjusted_cox),
    "\n")
cat("HR < 1:",
    sum(candidate_adjusted_cox$HR < 1),
    "\n")
cat("HR > 1:",
    sum(candidate_adjusted_cox$HR > 1),
    "\n")


# ============================================================
# 10. PROPORTIONAL HAZARDS DIAGNOSTIC
#     Top 20 adjusted Cox candidates
# ============================================================

top20_adjusted_genes <- head(
  candidate_adjusted_cox,
  20
)

ph_results <- vector(
  "list",
  nrow(top20_adjusted_genes)
)

for (i in seq_len(nrow(top20_adjusted_genes))) {

  gene <- top20_adjusted_genes$gene_symbol[i]

  gene_index <- which(
    gene_symbols == gene
  )[1]

  gene_expression <- as.numeric(
    expression_matrix[gene_index, ]
  )

  cox_data <- data.frame(
    OS_MONTHS = analysis_clinical$OS_MONTHS,
    OS_event = analysis_clinical$OS_event,
    AGE = analysis_clinical$AGE,
    SEX = analysis_clinical$SEX,
    expression = gene_expression
  )

  fit <- coxph(
    Surv(OS_MONTHS, OS_event) ~
      expression + AGE + SEX,
    data = cox_data
  )

  ph_test <- cox.zph(fit)

  expression_row <- which(
    rownames(ph_test$table) == "expression"
  )

  ph_results[[i]] <- data.frame(
    gene_symbol = gene,
    PH_p_value = ph_test$table[
      expression_row, "p"
    ]
  )
}

ph_results <- do.call(
  rbind,
  ph_results
)

ph_results$PH_assumption <- ifelse(
  ph_results$PH_p_value < 0.05,
  "Potential violation",
  "No evidence of violation"
)

ph_results <- ph_results[
  order(ph_results$PH_p_value),
]

write.csv(
  ph_results,
  file.path(
    results_dir,
    "top20_adjusted_cox_PH_diagnostic.csv"
  ),
  row.names = FALSE
)

cat("\nPH diagnostic:\n")
cat("Genes evaluated:", nrow(ph_results), "\n")
cat(
  "Potential violations:",
  sum(ph_results$PH_p_value < 0.05),
  "\n"
)

# ============================================================
# 11. Figure 1 — Top 20 survival-associated genes
#     Adjusted Cox regression
# ============================================================

top20_forest <- head(
  candidate_adjusted_cox,
  20
)

top20_forest$gene_symbol <- factor(
  top20_forest$gene_symbol,
  levels = rev(top20_forest$gene_symbol)
)

figure_1 <- ggplot(
  top20_forest,
  aes(
    x = HR,
    y = gene_symbol
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_errorbar(
    aes(
      xmin = CI_lower,
      xmax = CI_upper
    ),
    orientation = "y",
    width = 0.18,
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.8
  ) +
  scale_x_log10(
    breaks = c(0.2, 0.3, 0.5, 1, 2, 3, 5, 10),
    labels = c("0.2", "0.3", "0.5", "1", "2", "3", "5", "10")
  ) +
  labs(
    title = "Top 20 Survival-Associated Genes",
    subtitle = "Cox regression adjusted for age and sex",
    x = "Hazard Ratio (log scale)",
    y = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 12)
    ),
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 10)
    ),
    axis.text.x = element_text(
      size = 10
    ),
    axis.text.y = element_text(
      size = 10,
      face = "italic"
    ),
    axis.line = element_line(
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      linewidth = 0.5
    ),
    plot.margin = margin(
      12, 15, 12, 12
    )
  )

print(figure_1)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure_1_top20_adjusted_cox_forest_plot.png"
  ),
  plot = figure_1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
# ============================================================
# 12. FIGURE 2 — AGE–EXPRESSION ASSOCIATION
# ============================================================

age_results$neg_log10_FDR <- -log10(age_results$FDR)

age_results$significance <- ifelse(
  age_results$FDR < 0.05,
  "FDR < 0.05",
  "Not significant"
)

top_age_genes <- head(
  age_results[
    age_results$FDR < 0.05,
  ][order(
    age_results[
      age_results$FDR < 0.05,
    ]$FDR
  ), ],
  10
)

age_plot <- ggplot(
  age_results,
  aes(
    x = rho,
    y = neg_log10_FDR
  )
) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_point(
    aes(color = significance),
    size = 1.6
  ) +
  geom_text_repel(
    data = top_age_genes,
    aes(label = gene_symbol),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.7,
    point.padding = 0.35,
    min.segment.length = 0,
    force = 2,
    force_pull = 0.5
  ) +
  scale_color_manual(
    values = c(
      "FDR < 0.05" = "#2C7FB8",
      "Not significant" = "#BDBDBD"
    )
  ) +
  scale_shape_manual(
    values = c(
      "FDR < 0.05" = 16,
      "Not significant" = 1
    )
  ) +
  labs(
    title = "Gene Expression Associations with Age",
    subtitle = "Spearman correlation with age across 177 patients",
    x = "Spearman correlation (ρ)",
    y = expression(-log[10]("FDR")),
    color = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 12)
    ),
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 10)
    ),
    axis.text = element_text(
      size = 10
    ),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(
      size = 10
    ),
    plot.margin = margin(
      12, 15, 12, 12
    )
  )

print(age_plot)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure_2_age_expression_association.png"
  ),
  plot = age_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
# ============================================================
# 13. FIGURE 3 — SEX–EXPRESSION ASSOCIATION
# ============================================================

sex_results$neg_log10_FDR <- ifelse(
  sex_results$FDR == 0,
  35,
  -log10(sex_results$FDR)
)

sex_results$significance <- ifelse(
  sex_results$FDR < 0.05,
  "FDR < 0.05",
  "Not significant"
)

top_sex_genes <- head(
  sex_results[
    sex_results$FDR < 0.05,
  ][order(
    sex_results[
      sex_results$FDR < 0.05,
    ]$FDR
  ), ],
  10
)

sex_plot <- ggplot(
  sex_results,
  aes(
    x = effect_size,
    y = neg_log10_FDR
  )
) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_point(
    aes(color = significance),
    size = 1.6
  ) +
  geom_text_repel(
    data = top_sex_genes,
    aes(label = gene_symbol),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.7,
    point.padding = 0.35,
    min.segment.length = 0,
    force = 2,
    force_pull = 0.5
  ) +
  scale_color_manual(
    values = c(
      "FDR < 0.05" = "#2C7FB8",
      "Not significant" = "#BDBDBD"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 37),
    breaks = c(0, 5, 10, 15, 20, 25, 30, 35),
    expand = expansion(mult = c(0.02, 0.01))
  ) +
  labs(
    title = "Gene Expression Associations with Sex",
    subtitle = paste(
      "Wilcoxon rank-sum test across 177 patients",
      "FDR=0 values displayed at −log10(FDR)=35",
      sep = " · "
    ),
    x = "Rank-biserial effect size",
    y = expression(-log[10]("FDR")),
    color = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 12)
    ),
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 10)
    ),
    axis.text = element_text(
      size = 10
    ),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(
      size = 10
    ),
    plot.margin = margin(
      12, 15, 12, 12
    )
  )

print(sex_plot)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure_3_sex_expression_association.png"
  ),
  plot = sex_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# ============================================================
# 14. FIGURE 4 — OVERALL SURVIVAL KAPLAN–MEIER
# ============================================================

km_fit <- survfit(
  Surv(OS_MONTHS, OS_event) ~ 1,
  data = analysis_clinical
)

km_summary <- summary(
  km_fit,
  censored = TRUE
)

n_patients <- nrow(analysis_clinical)
n_events <- sum(analysis_clinical$OS_event)

median_os <- summary(km_fit)$table["median"]

km_data <- data.frame(
  time = km_summary$time,
  surv = km_summary$surv,
  lower = km_summary$lower,
  upper = km_summary$upper
)

km_plot <- ggplot(
  km_data,
  aes(
    x = time,
    y = surv
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15
  ) +
  geom_step(
    linewidth = 1,
    color = "#2C7FB8"
  ) +
  geom_point(
    data = data.frame(
      time = km_fit$time[km_fit$n.censor > 0],
      surv = km_fit$surv[km_fit$n.censor > 0]
    ),
    aes(x = time, y = surv),
    shape = 3,
    size = 2.2,
    color = "#2C7FB8"
  ) +
  geom_vline(
    xintercept = median_os,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = median_os + 1.5,
    y = 0.55,
    label = paste0(
      "Median OS = ",
      round(median_os, 2),
      " months"
    ),
    hjust = 0,
    size = 3.8
  ) +
  labs(
    title = "Overall Survival",
    subtitle = paste0(
      "Kaplan–Meier estimate · n = ",
      n_patients,
      " patients · ",
      n_events,
      " events"
    ),
    x = "Time (months)",
    y = "Overall survival probability"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 12)
    ),
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 10)
    ),
    axis.text = element_text(
      size = 10
    ),
    plot.margin = margin(
      12, 15, 12, 12
    )
  )

print(km_plot)

ggsave(
  filename = file.path(
    figures_dir,
    "Figure_4_overall_survival_KM.png"
  ),
  plot = km_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# ============================================================
# 15. FINAL SUMMARY
# ============================================================

cat("\n============================================\n")
cat("ANALYSIS COMPLETED\n")
cat("============================================\n")
cat("Genes:", n_genes, "\n")
cat("Samples:", n_samples, "\n")
cat("Age-associated genes (FDR < 0.05):",
    sum(age_results$FDR < 0.05, na.rm = TRUE),
    "\n")
cat("Sex-associated genes (FDR < 0.05):",
    sum(sex_results$FDR < 0.05, na.rm = TRUE),
    "\n")
cat("Adjusted Cox candidates:",
    nrow(candidate_adjusted_cox),
    "\n")
cat("PH violations among Top 20:",
    sum(ph_results$PH_p_value < 0.05),
    "\n")
cat("Median overall survival:",
    median_os,
    "months\n")
cat("============================================\n")
