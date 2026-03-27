source(file = "/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/R/00_functions.R")

# Counting geneset for expression and saving them in data/ora_archs4
counting_genesets(
  path_to_pd_objects = "data/paired_diff_archs4", 
  padj_score = "padj_expression", 
  out_dir = "data/ora_archs4/", 
  n_cores = 8
)

# Counting geneset for splicing and saving them in data/ora_archs4
counting_genesets(
  path_to_pd_objects = "data/paired_diff_archs4", 
  padj_score = "padj_expression", 
  out_dir = "data/ora_archs4/", 
  n_cores = 8
)