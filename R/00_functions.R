# Function to search metadata from ARCHS4 database h5-file - inputs described below
## *metadata_path* is the meta information such as title and characteristics columns from ARCHS4
## *keywords* are words chosen to be associated with the disease in question
## *max.distance* is the fraction specified as distance that each word that can differ from the original word (insertions, deletions, and substitutions).

meta_keyword_filter <- function(metadata, 
                                keywords, 
                                max.distance = 0.2) {
  
  # Keyword search in metadata
  metadata |> 
    dplyr::select(study, 
                  id, 
                  title, 
                  characteristics) |> 
    dplyr::filter(
      rowSums(sapply(keywords, function(keyword) {
        agrepl(pattern = keyword, 
               x = title, 
               ignore.case = TRUE, 
               max.distance = max.distance) |
          agrepl(pattern = keyword, 
                 x = characteristics, 
                 ignore.case = TRUE, 
                 max.distance = max.distance)
      })) > 0
    )
}

# Wrapper function to map all keywords chosen to search in metadata for each disease
data_mine_disease <- function(metadata_path,
                              keywords,
                              max.distance = 0.2,
                              disease_name = "disease_name") {
  
  # Reading metadata from path to filter studies where keywords are found
  metadata <- readr::read_rds(metadata_path)
  
  # index map to search for keywords in the metadata
  matching_studies_lists <- purrr::imap(
    keywords,
    ~ meta_keyword_filter(
      metadata = metadata,
      keywords = .x,
      max.distance = max.distance
    )
  )
  
  # Merging list rows to a dataframe
  matching_studies_df <- dplyr::bind_rows(matching_studies_lists)
  
  # Frequency table used to pull all relevant samples for each study in next step
  disease_freq_table <- table(matching_studies_df$study)
  
  # Additional extraction step that ensures all samples from extracted studies are pulled
  all_matching_metadata <- metadata |> 
    dplyr::filter(study %in% names(disease_freq_table)) |> 
    dplyr::select(study, id, title, characteristics) |> 
    dplyr::rename(gse_id = study,
                  gsm_id = id) |> 
    dplyr::mutate(condition = "")
  
  # write out the metadata that matches as excel file for manual curation
  writexl::write_xlsx(x = all_matching_metadata, path = paste0("data/metadata_archs4/", disease_name, "_metadata_filtered.xlsx"))
  
}


# Function to run pairedGSEA
## Function runs paired_diff for any kind of RNA-seq data with the new setup of having one comparisons file and one samples file.
## 
##
run_paired_diff <- function(n_cores = 4,
                            curated_comparison_file_path,
                            curated_samples_file_path,
                            curated_tcga_file_path = "data/curated_tcga/BRCA_LumA.xlsx",
                            archs4 = FALSE
                            ) 
{
  
  # Register cores
  if (archs4) {
    doMC::registerDoMC(cores = n_cores)
  }
    
  
  if (archs4){
    # Read in the two curated paired datasets "disease_comparisons" and "disease_samples"
    comparison <- readxl::read_xlsx(curated_comparison_file_path) |> 
      tibble::tibble() |> 
      dplyr::select(-notes) |> 
      dplyr::group_by(gse_id) |> 
      dplyr::mutate(row_id = dplyr::row_number()) |> 
      dplyr::ungroup()
    
    samples <- readxl::read_xlsx(curated_samples_file_path) |> 
      tibble::tibble() |> 
      dplyr::group_by(gse_id) |> 
      dplyr::mutate(row_id = dplyr::row_number()) |> 
      dplyr::ungroup()
    
  } else {
    samples_comparison <- readxl::read_xlsx(curated_tcga_file_path) |> 
      dplyr::rename(age = age_at_initial_pathologic_diagnosis,
             stage = ajcc_pathologic_tumor_stage) |> 
      dplyr::mutate(sample = gsub("-", "_", sample)) |> 
      tidyr::drop_na(age,
                     gender,
                     stage) |> 
      dplyr::filter(!stage %in% c("[Discrepancy]", "[Unknown]"))
    
      
  }
  
  if (archs4){
    # Joining the two curated paired datasets
    study_specific_merge <- dplyr::left_join(x = samples,
                                             y = comparison,
                                             by = c("gse_id", "row_id")) |> 
      dplyr::select(-row_id)
    
    # Curated samples 
    distinct_samples <- study_specific_merge |> 
      dplyr::pull(gsm_id)
    
    # Keeping a dataframe with only rows with gse_id, comparison, and Descriptions
    study_cond_desc_df <- study_specific_merge |> 
      tidyr::drop_na(comparison, 
                     description) |> 
      dplyr::select(gse_id,
                    comparison,
                    description)
    
  } else {
    # Curated samples 
    distinct_samples <- samples_comparison |> 
      dplyr::pull(sample)
    
    # Keeping a dataframe with only rows with gse_id, comparison, and Descriptions
    study_cond_desc_df <- samples_comparison |> 
      tidyr::drop_na(comparison, 
                     description) |> 
      dplyr::select(cancer_type_abbreviation,
                    comparison,
                    description)
  }
  
  # Collecting sample, gene, and transcript files
  ## samples 
  if (archs4) {
    archs4_samples_file <- "/home/databases/archs4/v2.latest/human_transcript_v2.latest.h5"
    samples <- rhdf5::h5read(file = archs4_samples_file,
                             name = "meta/samples/geo_accession")
  } else {
    if (!exists("tcgaIsoCount")) {
      load("/home/databases/tcga/tcga_Kallisto_est_counts.Rdata")
    }
    samples <- gsub("-", "_", colnames(tcgaIsoCount))
  }
  
  
  ## transcripts
  if (archs4){
    archs4_transcripts_file <- "/home/databases/archs4/v2.latest/human_tpm_v2.latest.h5"
    transcripts <- rhdf5::h5read(file = archs4_transcripts_file,
                                 name = "meta/transcripts/ensembl_id")
  } else {
    transcripts <- rownames(tcgaIsoCount)
  }
  
  ## genes
  if (archs4){
    archs4_genes_file <- "/home/databases/archs4/v2.latest/human_tpm_v2.latest.h5"
    genes <- rhdf5::h5read(file = archs4_genes_file,
                           name = "meta/transcripts/ensembl_gene")
  } else {
    
    genes <- readr::read_rds("/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/data/metadata_tcga/genes.rds")
  }
  
  # Filter and collect samples curated in the count matrices
  if (archs4){
    sample_locations <- which(samples %in% distinct_samples)
  } else {
    sample_locations <- which(samples %in% distinct_samples)
  }
  
  # Extracting expression data
  if (archs4){
    expression <- t(rhdf5::h5read(file = archs4_samples_file,
                                  name = "data/expression",
                                  index = list(sample_locations, 
                                               1:length(transcripts))))
  } else {
    expression <- tcgaIsoCount[, sample_locations, drop = FALSE] |> 
      dplyr::mutate(across(everything(), as.integer)) |> 
      as.matrix()
  }
  
  # Assigning rownames
  if (archs4) {
    exp_row_names <- character(nrow(genes))
  } else {
    exp_row_names <- character(length(transcripts))
  }
  
  
  # Iteration over rownames to combine gene and transcripts
  if (archs4){
    for (i in 1:nrow(genes)) {
      exp_row_names[i] <- paste0(genes[i],
                                 ":", 
                                 transcripts[i])
    }
  } else {
    for (i in 1:length(genes)) {
      exp_row_names[i] <- paste0(genes[i], ":", transcripts[i])
    }
  }
  
  # Assign the row names and column names to the data frame
  if (archs4){
    rownames(expression) <- exp_row_names
    colnames(expression) <- samples[sample_locations]
  } else {
    rownames(expression) <- exp_row_names
    colnames(expression) <- samples[sample_locations]
    colnames(expression) <- gsub("-", "_", colnames(expression))
      
  }
  
  # Subset expression matrix/dataframe 
  if (archs4){
    list_studies_merged_filt <- study_specific_merge |>
      dplyr::filter(gsm_id %in% colnames(expression))
  } else {
    list_studies_merged_filt <- samples_comparison |> 
      dplyr::filter(sample %in% colnames(expression)) 
  }
  
  # Paste gse_id and comparison together with separator ":"
  if (archs4){
    study_cond_desc_list <- split(
      study_cond_desc_df,
      stringr::str_c(study_cond_desc_df$gse_id, 
                     study_cond_desc_df$comparison, 
                     sep=' : ')
    )
  } else {
    study_cond_desc_list <- split(
      study_cond_desc_df,
      stringr::str_c(study_cond_desc_df$cancer_type_abbreviation, 
                     study_cond_desc_df$comparison, 
                     sep=' : ')
    )
  }
  # Running llply function to get a list from a list
  paired_diff_results <- plyr::llply(
    .data = study_cond_desc_list,
    .parallel = TRUE,
    .inform = TRUE,
    .progress = 'text',
    .fun = function(aDF) {
      
      ### For devel and debug
      # TCGA
      #aDF <- study_cond_desc_list[["BRCA : 6 vs 3"]]
      # ARCHS4
      #aDF <- study_cond_desc_list[["GSE116899 : 1 vs 2"]]
      
      
      ### Access sample information (incl subsetting)
      local_condition <- 
        stringr::str_split(aDF$comparison, ' vs ') |>
        unlist() |>
        as.integer()
      
      ### Finalizing subsets
      if (archs4){
        local_study_info <- list_studies_merged_filt |>
          dplyr::select(-comparison, 
                        -description) |>
          dplyr::filter(
            gse_id == aDF$gse_id,
            condition %in% local_condition
          ) 
      } else {
        local_study_info <- list_studies_merged_filt |>
          dplyr::select(-comparison, 
                        -description,
                        -Subtype_Immune_Model_Based,
                        -cancer_type_abbreviation) |>
          dplyr::filter(
            condition %in% local_condition) |> 
          dplyr::mutate(sample = gsub("-", "_", sample)) |> 
          tidyr::drop_na(age,
                         gender,
                         stage)
      }
      
      
      # Minimum count for each condition
      min_count <- min(table(local_study_info$condition))
      if(min_count < 2) {
        message("The minimum required of samples for each condition is 2, local_study_info do not have that.")
        return(NULL)
      }
      
      # Each subset/comparison needs to have 2 conditions check!
      cond_count <- length(unique(local_study_info$condition))
      if(cond_count != 2){
        message("There is not two different conditions in local_study_info")
        return(NULL)
      }
      
      # TCGA covariates flag
      covariate_cols <- c("age", "gender", "stage")
      
      covariate_arg <- if (all(covariate_cols %in% colnames(local_study_info))) {
        
        valid_covariates <- covariate_cols[
          sapply(covariate_cols, function(col) {
            n_unique <- length(unique(na.omit(local_study_info[[col]])))
            
            if (col == "gender") {
              n_unique > 1
            } else {
              n_unique > 5
            }
          })
        ]
        
        # Return NULL if none remain
        if (length(valid_covariates) > 0) {
          valid_covariates
        } else {
          message("No valid covariates (age, gender, stage) have more than one unique non-NA value.")
          NULL
        }
        
      } else {
        message("Not all required covariate columns (age, gender, stage) are present.")
        NULL
      }
      
      # Expression data archs4 or TCGA
      if (archs4){
        expression_obj <- expression[,local_study_info$gsm_id]
      } else {
        expression_obj <- expression[,local_study_info$sample]
      }
      
      # Sample column archs4 or TCGA
      if (archs4){
        samp_col <- 'gsm_id'
      } else {
        samp_col <- 'sample'
      }
      
      # logical used for running paired_diff on comparison with more than 20vs20.
      limma_flag <- nrow(local_study_info) > 30
      
      ### Running the paired_diff function from pairedGSEA
      paireDiff <- tryCatch({
        pairedGSEA::paired_diff(
          object = expression_obj, 
          group_col = 'condition', 
          sample_col = samp_col,
          baseline = local_condition[1],
          case = local_condition[2],
          metadata = local_study_info,
          use_limma = limma_flag,
          store_results = TRUE,
          experiment_title = aDF$description,
          covariates = covariate_arg
        )
      }, error = function(e) {
        message("paired_diff failed: ", conditionMessage(e))
        NULL
      })
      
      
      ### Return data
      return(paireDiff)
      
    }  
  )
}


# Function to run all tcga datasets with run_paired_diff
## 
## 
##
run_all_tcga <- function(tcga_dir, n_cores, out_dir) {
  files <- list.files(tcga_dir, pattern = "\\.xlsx$", full.names = TRUE)
  
  # Create directory only if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  for (f in files) {
    file_name <- tools::file_path_sans_ext(basename(f))
    out_file <- file.path(out_dir, paste0(file_name, "_paired_diff.rds"))
    
    # Skip if output file already exists
    if (file.exists(out_file)) {
      message("Skipping (already exists): ", basename(out_file))
      next
    }
    
    message("Running: ", basename(f))
    
    paired_diff_results <- run_paired_diff(
      n_cores = n_cores,
      curated_tcga_file_path = f,
      curated_samples_file_path = NULL,
      curated_comparison_file_path = NULL,
      archs4 = FALSE
    )
    
    readr::write_rds(paired_diff_results, file = out_file)
  }
}


# Function to run all archs4 datasets with run_paired_diff
## 
## 
##
run_all_archs4 <- function(datasets, base_dir, n_cores, out_dir) {
  
  # Create directory only if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  for (d in datasets) {
    
    sample_path <- file.path(base_dir, paste0(d, "_samples.xlsx"))
    comparison_path <- file.path(base_dir, paste0(d, "_comparisons.xlsx"))
    
    message("Running dataset: ", d)
    
    paired_diff_results <- run_paired_diff(
      n_cores = n_cores,
      curated_samples_file_path = sample_path,
      curated_comparison_file_path = comparison_path,
      curated_tcga_file_path = NULL,
      archs4 = TRUE
    )
    
    readr::write_rds(
      paired_diff_results,
      file = file.path(out_dir, paste0(d, "_paired_diff.rds"))
    )
  }
}


# Function convert list from pairedGSEA results to a tibble/dataframe
## 
## 
##

combine_lists_to_df <- function(pairedGSEA_list = paired_diff_results,
                                n_cores) {
  
  # Register n cores
  doMC::registerDoMC(cores = n_cores)
  
  # apply tibble class to lists
  paired_tibble <- lapply(pairedGSEA_list, tibble::as_tibble)
  
  # convert each list to dataframe
  plyr::ldply(
    .data = paired_tibble,
    .parallel = TRUE
  )
  
}


# Function count expression and splicing features
## 
## 
##
counting_features <- function(path = "data/paired_diff_archs4/", 
                              padj_score = "padj_splicing",
                              new_path = "data/counting_features_archs4/",
                              ...){
  
  # Create new directory if not present
  if (!dir.exists(new_path)) {
    dir.create(new_path)
  }
  
  # List all files to perform counting of features on
  files <- list.files(path = path, pattern = "\\.rds$", full.names = TRUE)
  
  # Walk across all files listed 
  purrr::walk(files, function(f) {
    
    # Converting lists to dataframe
    data <- readr::read_rds(f) |>
      combine_lists_to_df(...)
    
    # Counting features 
    result <- data |>
      dplyr::filter(.data[[padj_score]] < 0.05) |>
      dplyr::count(gene, name = "n") |>
      dplyr::arrange(desc(n))
    
    # Creating new filename
    new_name <- basename(f) |>
      stringr::str_remove("paired_diff") |>
      stringr::str_replace("\\.rds$", 
                           paste0("feature_counts_", padj_score, ".rds"))
    
    # New path with addition of new filename
    out_path <- file.path(new_path, new_name)
    
    # Write new RDS file
    readr::write_rds(result, out_path)
    
  })
}


# Function collects genesets and count both expression and splicing
## 
## 
##
counting_genesets <- function(path_to_pd_objects,
                              padj_score,
                              out_dir,
                              genesets = NULL,
                              ...) {
  
  # Register n cores
  doMC::registerDoMC(cores = n_cores)
  
  # Create output directory if needed
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  # List all input files
  files <- list.files(
    path = path_to_pd_objects, 
    pattern = "\\.rds$", 
    full.names = TRUE)
  
  # Read in genesets once
  if (is.null(genesets)) {
    genesets <- readr::read_rds(
      "/home/projects2/kvs_students/2026/sy_common_disease_iso/CrossDiseaseIsoformsR/data/genesets/genesets.rds"
    )
  }
  
  # Loop over files
  purrr::walk(files, function(f) {
    
    # Parsing old name to set new name
    new_base_fname <- gsub(
      "paired_diff",
      "",
      tools::file_path_sans_ext(basename(f))
    )
    
    # Read in data and combine to dataframe
    data <- readr::read_rds(file = f) |>
      combine_lists_to_df()
    
    pairedOraResults <- plyr::dlply(
      .data = data,
      .variables = c(".id"),
      .parallel = TRUE,
      .fun = function(data) {
        pairedGSEA::paired_ora(
          paired_diff_result = data,
          gene_sets = genesets
        )
      }
    ) |>
      combine_lists_to_df(...)
    
    ora_results <- pairedOraResults |>
      dplyr::filter(.data[[padj_score]] < 0.05) |>
      dplyr::count(pathway, name = "n") |>
      dplyr::arrange(desc(n))
    
    # Build output file path
    out_file <- file.path(out_dir, paste0(new_base_fname, "ora_results.rds"))
    
    # Save result
    saveRDS(ora_results, file = out_file)
    
    message("Saved: ", out_file)
  })
}

