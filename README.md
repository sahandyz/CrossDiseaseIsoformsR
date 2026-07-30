# CrossDiseaseIsoformsR
This repository contains the R scripts used throughout my thesis, from the paired differential analysis, meta-analyses and result simplification, and drug target discovery.

## Workflow Overview

The scripts have been setup to be used in chronological order to replicate all the analyses of this thesis.

| Script | Description |
|--------|-------------|
| **00_functions.R** | Contains all custom functions used throughout the thesis. Source this script before running the analysis. |
| **01_run_archs4_pairedGSEA.R** | Performs paired differential expression and splicing analyses on the ARCHS4 datasets. |
| **02_run_tcga_pairedGSEA.R** | Performs paired differential expression and splicing analyses on the TCGA datasets. |
| **03_prepare_genesets.R** | Prepares the gene-sets used for both Over-Representation Analysis (ORA) and Functional Class Scoring (FCS) throughout the thesis. |
| **04_preliminary_visualizations.R** | Generates exploratory visualizations, including boxplots illustrating the importance of isoforms compared to genes. |
| **05_cutoff.R** | Determines frequency-based thresholds for genes and isoforms. Features were required to appear in up to **9 diseases** and in at least **5–6 studies** within each disease. |
| **06_fgsea_cutoff.R** | Performs Over-Representation Analysis (ORA) using `fgsea::fora()`. |
| **07_geneSetSimplifyR_cutoff.R** | Creates `geneSetSimplifyR` objects from the ORA results and generates visualizations to simplify interpretation of high dimensional gene-set results. |
| **08_gene_level_analysis.R** | Produces gene-level volcano plots and isoform-level splicing scatter plots, highlighting disease-associated genes and genes that spliced. |
| **09_formal_p_value_integration.R** | Utilizes statistical evidence by applying **Fisher's Method** across studies for each disease, followed by **Edgington's Method** across diseases. Performs Functional Class Scoring (FCS) using `fgseaMultilevel()` and prepares `geneSetSimplifyR` objects for downstream visualization. |
| **10_geneSetSimplifyR_p_val_int.R** | Uses the `geneSetSimplifyR` objects generated in Script 09 to visualize and simplify the FCS enrichment results. |
| **11_diseases_JensenLab.R** | Performs FCS using the JensenLab **DISEASES** gene-set collection for all diseases except **Tuberculosis**, as the corresponding genes were unavailable. |
| **12_Cross_Indication_Target_Discovery.R** | Performs drug repurposing by combining **leadingEdge genes** from FCS and **overlapGenes** from ORA. These gene-sets are intersected against the **ChEMBL human target database**, and evidence from both analyses are integrated to identify high-confidence therapeutic targets with known manually curated bioactive molecules. |

## Analysis Summary

The overall workflow consists of:

1. Differential expression and splicing analyses.
2. Gene-set preparation.
3. Exploratory visualization of Genes vs Isoforms.
4. Frequency-based feature filtering and summarization.
5. Over-Representation Analysis (ORA).
6. Simplification and visualization of the ORA using `geneSetSimplifyR`.
7. Cross-study and cross-disease p-value integration.
8. Functional Class Scoring (FCS).
9. Simplification and visualization of the FCS using `geneSetSimplifyR`.
10. Gene- and isoform-level visualization.
11. Disease-specific enrichment using JensenLab DISEASES.
12. Cross-indication drug target discovery using ChEMBL.

## Main R Packages

The analysis relies on several R packages, including:

- **DESeq2**
- **DEXSeq**
- **fgsea**
- **geneSetSimplifyR**
- **msigdbr**

## R Packages

The analysis was performed using the following R packages (some of which use the **Main R Packages**):

| Package | Version |
|:--------|:-------:|
| data.table | 1.16.0 |
| doMC | 1.3.8 |
| dplyr | 1.1.4 |
| fgsea | 1.30.0 |
| forcats | 1.0.0 |
| ggbeeswarm | 0.7.3 |
| ggplot2 | 3.5.1 |
| ggrepel | 0.9.6 |
| ggridges | 0.5.6 |
| httr | 1.4.7 |
| jsonlite | 1.8.8 |
| metap | 1.14 |
| pairedGSEA | 1.13.0 |
| patchwork | 1.2.0 |
| plyr | 1.8.9 |
| purrr | 1.0.2 |
| readr | 2.1.5 |
| readxl | 1.4.3 |
| rhdf5 | 2.48.0 |
| rlang | 1.1.4 |
| stats | 4.4.1 |
| stringr | 1.5.1 |
| tibble | 3.2.1 |
| tidyr | 1.3.1 |
| tools | 4.4.1 |
| viridisLite | 0.4.2 |
| writexl | 1.5.0 |

These packages were used throughout the analysis. Additional package dependencies may be loaded within individual scripts where required.