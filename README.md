# Ovarian cancer integrated miRNA-mRNA analysis code

This repository contains analysis scripts and reproducibility materials for the manuscript:

Integrated miRNA-mRNA network analysis identifies a miR-15a-5p-centered angiogenic program and biomarker-linked regulatory modules in ovarian cancer.

## Contents

- scripts/00_manuscript_readiness_check.R  
  Checks the presence of expected manuscript tables, figures, supplementary files, and analysis outputs.

- scripts/06_TCGA_module5_covariate_sensitivity.R  
  Exploratory TCGA-OV sensitivity analysis for the Module 5 miR-15a-5p angiogenesis association.

- scripts/legacy_miRNA_mRNA_integration_review_before_use.R  
  Earlier miRNA-mRNA integration script retained for provenance review. This file should not be treated as the final revised network-construction script unless reconciled with the final exported network table.

- sessionInfo.txt  
  R session information from the revision environment.

## Important note

The final manuscript network is described as an expression-inferred network. The final exported edge table contains 4,614 inferred miRNA-mRNA edges linking 24 miRNAs and 1,236 mRNAs. All retained final edges are annotated as novel_inferred, and the final edge score is defined in the revised Methods.

## Data availability

Raw sequencing and expression data are deposited separately. Large raw files such as FASTQ, RCC, BAM, and full expression matrices are not included in this code repository.
