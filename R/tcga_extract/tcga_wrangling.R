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
survival_dat <- readr::read_delim("/home/databases/tcga/Survival_SupplementalTable_S1_20171025_xena_sp.txt") |> 
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
  left_join(immune_substype, by = "sample") |> 
  drop_na(age_at_initial_pathologic_diagnosis,
          gender,
          ajcc_pathologic_tumor_stage)

# Check above dataframe
check_relevant_sample <- relevant_cancer_samples |> 
  filter(subtype_selected %in% c("BRCA.Basal",
                                 "BRCA.Her2",
                                 "BRCA.LumA",
                                 "BRCA.LumB",
                                 "BRCA.Normal",
                                 "GI.CIN",
                                 "GI.GS",
                                 "GI.HM-indel",
                                 "HNSC.Atypical",
                                 "HNSC.Basal",
                                 "HNSC.Classical",
                                 "HNSC.Mesenchymal",
                                 "LUSC.basal",
                                 "LUSC.classical",
                                 "LUSC.primitive",
                                 "LUSC.secretory",
                                 "GI.CIN",
                                 "GI.EBV",
                                 "GI.GS",
                                 "GI.HM-indel"),
         !sample_type_id == "11")


get_cancers_no_sub <- relevant_cancer_samples |> 
  filter(subtype_selected %in% c(NA, 
                                 "KICH.Eosin.0",
                                 "KIRC.1",
                                 "KIRC.2",
                                 "KIRC.3",
                                 "KIRC.4",
                                 "KIRC.NA",
                                 "KIRP.C1",
                                 "KIRP.C2a",
                                 "LIHC.iCluster:1",
                                 "LIHC.iCluster:2",
                                 "LIHC.iCluster:3",
                                 "LUAD.2",
                                 "LUAD.3",
                                 "LUAD.4",
                                 "LUAD.5",
                                 "LUAD.6",
                                 "THCA.1",
                                 "THCA.2",
                                 "THCA.3",
                                 "THCA.4",
                                 "THCA.5"),
         !sample_type_id == "11") |> 
  drop_na(Subtype_Immune_Model_Based) |> 
  mutate(subtype_selected = case_when(
    Subtype_Immune_Model_Based == "Wound Healing (Immune C1)" ~ paste0(cancer_type_abbreviation, ".C1"),
    Subtype_Immune_Model_Based == "IFN-gamma Dominant (Immune C2)" ~ paste0(cancer_type_abbreviation, ".C2"),
    Subtype_Immune_Model_Based == "Inflammatory (Immune C3)" ~ paste0(cancer_type_abbreviation, ".C3"),
    Subtype_Immune_Model_Based == "Lymphocyte Depleted (Immune C4)" ~ paste0(cancer_type_abbreviation, ".C4"),
    Subtype_Immune_Model_Based == "Immunologically Quiet (Immune C5)" ~ paste0(cancer_type_abbreviation, ".C5"),
    Subtype_Immune_Model_Based == "TGF-beta Dominant (Immune C6)" ~ paste0(cancer_type_abbreviation, ".C6"),
  ))

healthy_samples <- relevant_cancer_samples |> 
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
  ) & sample_type_id == 11)

merge_cancers <- bind_rows(get_cancers_no_sub, 
                           check_relevant_sample,
                           healthy_samples)

final_cancer_df <- merge_cancers |> 
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
         ),
         subtype_selected = case_when(
           subtype_selected == "GI.CIN" ~ paste0(cancer_type_abbreviation, "_", subtype_selected),
           subtype_selected == "GI.GS" ~ paste0(cancer_type_abbreviation, "_", subtype_selected),
           subtype_selected == "GI.HM-indel" ~ paste0(cancer_type_abbreviation, "_", subtype_selected),
           subtype_selected == "GI.EBV" ~ paste0(cancer_type_abbreviation, "_", subtype_selected),
           TRUE ~ subtype_selected
         ),
         subtype_selected = str_replace(subtype_selected, "\\.", "_")) |>
  arrange(condition) |> 
  drop_na(age_at_initial_pathologic_diagnosis,
          gender,
           ajcc_pathologic_tumor_stage) |>
  group_by(cancer_type_abbreviation,
           sample_type_id,
           subtype_selected) |>
  reframe(n = n()) |>
  filter(n > 24)


# write out excel for final_cancer_df
writexl::write_xlsx(x = final_cancer_df,
                    path = "data/metadata_tcga/all_cancer_subtypes_11MAY.xlsx")

# Overview
x <- final_cancer_df |> 
  group_by(cancer_type_abbreviation,
           subtype_selected,
           sample_type_id,
           condition) |>
  reframe(n = n()) |> 
  arrange(condition)

# Get correct gene and isoform annotations for TCGA data.
# Loading annotation of transcript and genes
load("/home/databases/tcga/gencode.v23.chr_patch_hapl_scaff.annotation.Rdata")

# Creating vectorized Dframe
anno <- S4Vectors::mcols(gen23)

# Pulling only ENST and ENSG IDs from annotation
enst_ensg_annotation <- unique(data.frame(
  transcript_id = anno$transcript_id,
  gene_id = anno$gene_id,
  stringsAsFactors = FALSE
))

# Dropping NA transcripts
enst_ensg_filtered <- enst_ensg_annotation[!is.na(enst_ensg_annotation$transcript_id), ]

# Matching gene IDs and transcript IDs
matched_gene_id <- enst_ensg_filtered$gene_id[match(transcripts, enst_ensg_filtered$transcript_id)]

saveRDS(matched_gene_id, "data/metadata_tcga/genes.rds")
