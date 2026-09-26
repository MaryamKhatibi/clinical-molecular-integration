# Clinical–Molecular Integration in Pancreatic Adenocarcinoma

## Overview
This project integrates transcriptomic gene-expression data with clinical characteristics and overall-survival information in pancreatic adenocarcinoma.

**Source:** cBioPortal — Pancreatic Adenocarcinoma (TCGA, PanCancer Atlas)  
**Study ID:** `paad_tcga_pan_can_atlas_2018`

### Research question
> Is gene expression associated with clinical characteristics and patient survival in a pancreatic adenocarcinoma cohort?

## Cohort
- Integrated samples: 177
- Genes/features after preprocessing: 19,127
- Overall-survival events: 92
- Censored observations: 85

Seven clinical samples lacked matching expression samples and were excluded from the integrated analysis.

## Analysis workflow
1. Clinical–molecular data integration
2. Expression quality control
3. Log2(x + 1) transformation
4. Low-detection gene filtering
5. Spearman association with age
6. Wilcoxon rank-sum association with sex
7. Kaplan–Meier overall-survival analysis
8. Gene-wise unadjusted Cox regression
9. Gene-wise Cox regression adjusted for age and sex
10. Benjamini–Hochberg FDR correction
11. Proportional-hazards diagnostics
12. Reproducible figure generation

## Main figures
- Figure 1 — Top 20 survival-associated genes after adjustment for age and sex
- Figure 2 — Gene-expression associations with age
- Figure 3 — Gene-expression associations with sex
- Figure 4 — Overall Kaplan–Meier survival curve

## Methodological notes
Expression values were treated as continuous RSEM-like measurements. Downstream association and survival analyses used log2(x + 1)-transformed expression values.

Adjusted Cox hazard ratios describe associations with overall survival after adjustment for age and sex. They are observational associations and should not be interpreted as causal effects.

A stage-adjusted Cox analysis was explored but was not retained as a main analysis because the available stage distribution was highly imbalanced.

## Repository structure
```text
clinical-molecular-integration/
├── data/
│   └── README.md
├── scripts/
│   ├── 01_data_preparation.R
│   ├── 02_expression_processing.R
│   ├── 03_build_analysis_dataset.R
│   ├── 04_expression_clinical_integration.R
│   ├── 05_expression_preprocessing.R
│   └── 06_expression_clinical_association.R
├── figures/
├── results/
│   ├── tables/
│   └── preprocessed/
├── docs/
├── .gitignore
└── README.md
```

## Reproducibility
The final analysis is implemented in `scripts/06_expression_clinical_association.R`. Earlier scripts prepare and validate the inputs. The final script performs the statistical analyses and generates the four main figures.

## Data source
cBioPortal: https://www.cbioportal.org/

Study: **Pancreatic Adenocarcinoma (TCGA, PanCancer Atlas)**

Raw source data and large intermediate objects are not included in this repository.

## Software
- R
- survival
- ggplot2
- ggrepel

## Disclaimer
This repository contains an independent computational reanalysis of publicly available molecular and clinical data. Findings are exploratory associations and are not intended to establish clinical utility or causality.
