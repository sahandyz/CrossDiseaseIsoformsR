# Argument from SLURM
args <- commandArgs(trailingOnly = TRUE)

# Dataset index argument
dataset <- args[1]

source("/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/R/00_functions.R")

# Running all archs4 diseases.
run_all_archs4(
  datasets = c(dataset),
  base_dir = "/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/data/curated_archs4/",
  n_cores = 10,
  out_dir = "/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/data/pairedV2_diff_archs4"
)
