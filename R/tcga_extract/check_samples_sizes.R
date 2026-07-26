library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)



# KICH healthy goes from 25 to 18 - REMOVE
# KIRC_C4 goes from 27 to 23 - REMOVE
# LUAD_C6 goes from 27 to 16 - REMOVE
# LUSC_Primitive goes from 26 to 16 - REMOVE
# STAD_GI_EBV goes from 30 to 14 - REMOVE



filter_one_file <- function(file) {
  dat <- readxl::read_excel("data/curated_tcga_immune/THCA_C3.xlsx")
  unique(dat$condition)
  local_condition <- c("49")  # adjust this
  
  before_n <- nrow(dat)
  
  filtered <- dat |>
    dplyr::select(
      -comparison,
      -description,
      -Subtype_Immune_Model_Based,
      -cancer_type_abbreviation
    ) |>
    dplyr::filter(condition %in% local_condition) |>
    dplyr::mutate(sample = gsub("-", "_", sample)) |>
    tidyr::drop_na(age_at_initial_pathologic_diagnosis, gender, ajcc_pathologic_tumor_stage) |>
    dplyr::add_count(gender, name = "n_gender") |>
    dplyr::add_count(ajcc_pathologic_tumor_stage, name = "n_stage") |>
    dplyr::filter(n_gender >= 5, n_stage >= 5) |>
    dplyr::select(-n_gender, -n_stage)
  
  tibble(
    file = basename(file),
    samples_before = before_n,
    samples_after = nrow(filtered),
    samples_removed = before_n - nrow(filtered)
  )
}

results <- list.files(
  folder_path,
  pattern = "\\.xlsx?$",
  full.names = TRUE
) |>
  purrr::map_dfr(filter_one_file)

results
