# Sourcing functions
source("R/00_functions.R")

library(data.table)
library(dplyr)
library(fgsea)

# url for jensenlab
url <- "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv"

dis <- fread(url, header = FALSE)

# First 4 columns are always:
# gene identifier, gene name, disease identifier, disease name was needed
colnames(dis)[1:4] <- c("gene_id", "gene_symbol", "disease_id", "disease_name")

# inspect disease names to see what there and whats not
grep("arthritis|alzheimer's disease|2 diabete|Hypertension|Familial atrial fibrillation|Cardiomyopathy|Dilated cardiomyopathy|Restrictive cardiomyopathy|Aortic aneurysm|Thoracic aortic aneurysm|Thrombophilia|Essential thrombocythemia|Polyarteritis nodosa|Bronchiectasis|crohn's disease|ulcerative colitis|inflammatory bowel disease|obesity|tuberculosis|cancer|Mycobacterium tuberculosis", unique(dis$disease_name),
     ignore.case = TRUE, value = TRUE)


cvd_diseases <- c(
  "Hypertension",
  "Familial atrial fibrillation",
  "Cardiomyopathy",
  "Dilated cardiomyopathy",
  "Restrictive cardiomyopathy",
  "Aortic aneurysm",
  "Thoracic aortic aneurysm",
  "Thrombophilia",
  "Essential thrombocythemia",
  "Polyarteritis nodosa"
)

# Creating the list of gene-sets
jensen_pathways <- list(
  AD = dis |>
    filter(tolower(disease_name) == "alzheimer's disease") |>
    pull(gene_symbol) |>
    unique(),
  
  COPD = dis |>
    filter(tolower(disease_name) == "bronchiectasis") |>
    pull(gene_symbol) |>
    unique(),
  
  CVD = dis |>
    filter(tolower(disease_name) %in% tolower(cvd_diseases)) |>
    pull(gene_symbol) |>
    unique(),
  
  IBD = dis |>
    filter(tolower(disease_name) %in% c("crohn's disease",
                                        "ulcerative colitis",
                                        "inflammatory bowel disease")) |>
    pull(gene_symbol) |>
    unique(),
  
  OB = dis |>
    filter(tolower(disease_name) == "obesity") |>
    pull(gene_symbol) |>
    unique(),
  
  RA = dis |>
    filter(tolower(disease_name) == "rheumatoid arthritis") |>
    pull(gene_symbol) |>
    unique(),
  
  T2D = dis |>
    filter(tolower(disease_name) %in% c("type 2 diabetes mellitus",
                                        "type-2 diabetes")) |>
    pull(gene_symbol) |>
    unique(),
  
  CANCERS = dis |> 
    filter(tolower(disease_name) %in% c("cancer",
                                        "urinary bladder cancer",
                                        "breast cancer",
                                        "lung cancer")) |> 
    pull(gene_symbol) |>
    unique()
)

# new disease list with ensg-id to hgmc translation
new_dfs <- readRDS("data/ensg_to_symbol/new_dfs.rds")

# aggregation p-values
meta_per_disease <- meta_pvalues_by_gene(new_dfs)

# ranking of genes for FCS
make_ranks <- function(df, p_col) {
  df |>
    dplyr::mutate(
      p = pmax(.data[[p_col]], .Machine$double.xmin),
      rank = -log10(p)
    ) |>
    dplyr::select(gene, rank) |>
    tibble::deframe() |>
    sort(decreasing = TRUE)
}

# Applying ranks to all diseases
expr_fisher_ranks <- lapply(
  meta_per_disease,
  make_ranks,
  p_col = "expression_meta_p"
)

# Fgsea function to run with multiLevel
run_fgsea_matched <- function(ranks_list, pathways, minsize = 1) {
  common_diseases <- intersect(names(ranks_list), names(pathways))
  
  purrr::map_dfr(common_diseases, function(dis) {
    fgsea::fgseaMultilevel(
      pathways = pathways[dis],
      stats = ranks_list[[dis]],
      minSize = minsize
    ) |>
      dplyr::arrange(padj) |>
      dplyr::mutate(disease = dis)
  })
}

# Running multilevel on jensenLab Gene-sets and getting a data frame back with all
expr_fgsea_results <- run_fgsea_matched(
  ranks_list = expr_fisher_ranks,
  pathways = jensen_pathways
)


AD_expr_jensen <- expr_fgsea_results$AD
COPD_expr_jensen <- expr_fgsea_results$COPD
CVD_expr_jensen <- expr_fgsea_results$CVD
IBD_expr_jensen <- expr_fgsea_results$IBD
OB_expr_jensen <- expr_fgsea_results$OB
RA_expr_jensen <- expr_fgsea_results$RA
TB_expr_jensen <- expr_fgsea_results$TB
T2D_expr_jensen <- expr_fgsea_results$T2D
CANCERS_expr_jensen <- expr_fgsea_results$CANCERS

sp_fgsea_results <- run_fgsea_multilevel_dis(
  ranks_list = sp_fisher_ranks,
  pathways = jensen_pathways
)

AD_sp_jensen <- sp_fgsea_results$AD
