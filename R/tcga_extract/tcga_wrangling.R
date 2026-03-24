library(readr)
library(dplyr)
library(janitor)
library(tidyr)
library(stringr)

# Load isoform counts - tcgaIsoCount
#load("/home/databases/tcga/tcga_Kallisto_est_counts.Rdata") 

# Load isoform and gene info
#load("/home/databases/tcga/gencode.v23.chr_patch_hapl_scaff.annotation.Rdata")

# Load Clusters
clusters <- read_delim("/home/databases/tcga/TCGA_PanCan33_iCluster_k28_tumor.gz") |> 
  janitor::clean_names()

# Load phenotypes
pheno <- read_tsv("/home/databases/tcga/TCGA_phenotype_denseDataOnlyDownload.tsv.gz") |> 
  janitor::clean_names()

# Load subtypes
subtypes <- read_tsv("/home/databases/tcga/TCGASubtype.20170308.tsv.gz") |> 
  janitor::clean_names()

# Load survival data
survival_dat <- read_delim("/home/databases/tcga/Survival_SupplementalTable_S1_20171025_xena_sp.txt") |> 
  janitor::clean_names()

# Load Immune subtype file
immune_substype <- read_delim("/home/databases/tcga/Subtype_Immune_Model_Based.txt.gz")

# subset cols for pheno, subtypes, and survival
subset_pheno <- pheno |> 
  select(sample, 
         sample_type_id)

subset_subtypes <- subtypes |> 
  select(sample_id, 
         subtype_selected)

subset_survival <- survival_dat |> 
  select(sample, 
         cancer_type_abbreviation, 
         age_at_initial_pathologic_diagnosis, 
         gender, 
         ajcc_pathologic_tumor_stage) |> 
  filter(cancer_type_abbreviation %in% c("BRCA",
                                         "COAD",
                                         "HNSC",
                                         "KICH",
                                         "KIRC",
                                         "KIRP",
                                         "LIHC",
                                         "LUAD",
                                         "LUSC",
                                         "PRAD",
                                         "STAD",
                                         "THCA",
                                         "UCEC"
  ))

# join the three subset dataframes
relevant_cancer_samples <- subset_survival |> 
  left_join(subset_subtypes, by = c("sample" = "sample_id")) |> 
  left_join(subset_pheno, by = "sample") |> 
  left_join(immune_substype, by = "sample")

# Check above dataframe
check_relevant_sample <- relevant_cancer_samples |> 
  group_by(cancer_type_abbreviation,
           sample_type_id,
           subtype_selected) |> 
  reframe(n = n()) |> 
  filter(n > 24) 



final_cancer_df <- relevant_cancer_samples |> 
  mutate(subtype_selected = str_replace(subtype_selected,
                                        "\\.(\\d+)",
                                        ".C\\1"),
         subtype_selected = str_replace(subtype_selected,
                                        "\\.iCluster:(\\d+)",
                                        ".C\\1")) |> 
  group_by(cancer_type_abbreviation,
           sample_type_id,
           subtype_selected) |> 
  add_count(name = "n") |>  
  filter(n > 24,
         !grepl("\\.NA", subtype_selected),
         !(is.na(subtype_selected) & sample_type_id == "01")) |> 
  select(-n) |> 
  ungroup() |> 
  group_by(cancer_type_abbreviation,
           subtype_selected,
           sample_type_id) |> 
  mutate(condition = cur_group_id(),
         subtype_selected = case_when(
           (is.na(subtype_selected) | subtype_selected == "BRCA.Normal") & sample_type_id == "11" ~ paste0(cancer_type_abbreviation, "_Control"),
           TRUE ~ subtype_selected
         )) |>
  arrange(condition)

# Overview
x <- final_cancer_df |> 
  group_by(cancer_type_abbreviation,
           subtype_selected,
           sample_type_id,
           condition) |>
  reframe(n = n()) |> 
  arrange(condition)

# write out excel for final_cancer_df
writexl::write_xlsx(x = final_cancer_df,
                    path = "data/metadata_tcga/all_cancer_subtypes.xlsx")