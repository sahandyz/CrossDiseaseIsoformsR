library(readr)
library(dplyr)
library(janitor)

# Load isoform counts
#load("/home/databases/tcga/tcga_Kallisto_est_counts.Rdata")

# Load isoform and gene info
#load("/home/databases/tcga/gencode.v23.chr_patch_hapl_scaff.annotation.Rdata")

# Load phenotypes
pheno <- read_tsv("/home/databases/tcga/TCGA_phenotype_denseDataOnlyDownload.tsv.gz") %>% 
  janitor::clean_names()

# Load subtypes
subtypes <- read_tsv("/home/databases/tcga/TCGASubtype.20170308.tsv.gz") %>% 
  janitor::clean_names()

subtypes_subset <- subtypes
  

# Load survival data
survival_dat <- read_delim("/home/databases/tcga/Survival_SupplementalTable_S1_20171025_xena_sp.txt") %>% 
  janitor::clean_names()

# Keep relevant cols in survival data
survival_dat_subset <- survival_dat %>% 
  select(sample, 
         cancer_type_abbreviation, 
         age_at_initial_pathologic_diagnosis, 
         gender, 
         ajcc_pathologic_tumor_stage)

# Extract cancers with more than 24 healthy controls
above25_pheno <- pheno %>% 
  group_by(sample_type_id, primary_disease) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  filter(sample_type_id == "11",
         n > 24) 

# Subset pheno types of 12 cancer types with >24 healthy control (Based on Krivi study)
subset_pheno <- pheno %>% 
  filter(
    primary_disease %in% c(
      "breast invasive carcinoma",
      "colon adenocarcinoma",
      "head & neck squamous cell carcinoma",
      "kidney chromophobe",
      "kidney clear cell carcinoma",
      "kidney papillary cell carcinoma",
      "liver hepatocellular carcinoma",
      "lung adenocarcinoma",
      "lung squamous cell carcinoma",
      "prostate adenocarcinoma",
      "stomach adenocarcinoma",
      "thyroid carcinoma",
      "uterine corpus endometrioid carcinoma"
    )
  )

# Keep all sample from subset_pheno which are also in survival_dat_subset
good_enough_cancers <- inner_join(subset_pheno, 
                                  survival_dat_subset, 
                                  by = "sample")

# Healthy samples from good enoguh 
healthy_good_cancers <- good_enough_cancers %>% 
  filter(sample_type_id == "11")

# Disease samples from good enough 
disease_good_cancers <- good_enough_cancers %>% 
  filter(sample_type_id != "11")

# Adding the subtypes on good_enough_cancers
good_cancers_subtypes <- left_join(good_enough_cancers, 
                                    subtypes_subset, 
                                    by = c("sample" = "sample_id"))

# Check to see how many subtypes have 25 samples
good_enough_subtypes <- good_cancers_subtypes %>% 
  group_by(subtype_m_rna, cancer_type_abbreviation) %>% 
  summarize(n = n()) %>% 
  arrange(cancer_type_abbreviation) %>% 
  filter()
  


# Dieases specific data saving with subtypes
