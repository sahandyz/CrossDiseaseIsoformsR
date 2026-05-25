source("R/00_functions.R")

# Counting features for expression and saving them in data/counting_features_archs4
counting_features(
  path = "data/paired_diff_archs4/",
  padj_score = "padj_expression", 
  new_path = "data/counting_features_archs4",
  n_cores = 8
)

# Counting features for splicing and saving them in data/counting_features_archs4
counting_features(
  path = "data/paired_diff_archs4/",
  padj_score = "padj_splicing", 
  new_path = "data/counting_features_archs4",
  n_cores = 8
)
