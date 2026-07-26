# ---- Preliminary RNA-seq analysis functions ----

# Function to search metadata from ARCHS4 database h5-file - inputs described below
## *metadata_path* is the meta information such as title and characteristics columns from ARCHS4
## *keywords* are words chosen to be associated with the disease in question
## *max.distance* is the fraction specified as distance that each word that can differ from the original word (insertions, deletions, and substitutions).

meta_keyword_filter <- function(
    metadata, 
    keywords, 
    max.distance = 0.2
) 
{
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
data_mine_disease <- function(
    metadata_path,
    keywords,
    max.distance = 0.2,
    disease_name = "disease_name"
) 
{
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
    dplyr::select(study, 
                  id, 
                  title, 
                  characteristics) |> 
    dplyr::rename(gse_id = study,
                  gsm_id = id) |> 
    dplyr::mutate(condition = "")
  
  # write out the metadata that matches as excel file for manual curation
  writexl::write_xlsx(x = all_matching_metadata, 
                      path = paste0("data/metadata_archs4/", 
                                    disease_name, "_metadata_filtered.xlsx"))
  
}

# filter_low_levels_of_covariates <- function(metadata, var, group_col) {
#   tab <- table(metadata[[var]], metadata[[group_col]])
#   
#   # Levels where any condition has ≤ 1 sample
#   bad_levels <- rownames(tab)[apply(tab, 1, function(x) any(x <= 1))]
#   
#   metadata[!metadata[[var]] %in% bad_levels, , drop = FALSE]
# }

# Function to run pairedGSEA
## Function runs paired_diff for ARCHS4 or TCGA RNA-seq data with - the setup being: having one comparisons file and one samples file for ARCHS4.
## For the TCGA part of the code the columns needed are: gender, stage,	subtype_selected, sample_type_id, Subtype_Immune_Model_Based, condition, comparison, and description.
## number of cores, relevant paths, and whether the data needed to be run is ARCHS4 or TCGA can be specified with archs4 = TRUE (for ARCHS4), and archs4 = FALSE (for TCGA)
## The data is pulled from relevant files, prepared and run in parallel if needed with paired_diff from pairedGSEA, with multiple checks and filters.
run_paired_diff <- function(
    n_cores = 12,
    curated_comparison_file_path,
    curated_samples_file_path,
    curated_tcga_file_path = "data/curated_tcga/LIHC_C3.xlsx",
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
    samples_comparison <- readxl::read_xlsx(curated_tcga_file_path) 
    
    # Keeping a dataframe with only rows with gse_id, comparison, and Descriptions
    study_cond_desc_df <- samples_comparison |> 
      tidyr::drop_na(comparison, 
                     description) |> 
      dplyr::select(cancer_type_abbreviation,
                    comparison,
                    description)
    
    samples_comparison <- samples_comparison |> 
      dplyr::rename(age = age_at_initial_pathologic_diagnosis,
                    stage = ajcc_pathologic_tumor_stage) |> 
      dplyr::mutate(sample = gsub("-", 
                                  "_", 
                                  sample)) |> 
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
      dplyr::mutate(across(everything(), 
                           as.integer)) |> 
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
      #aDF <- study_cond_desc_list[["LIHC : 29 vs 28"]]
      # ARCHS4
      #aDF <- study_cond_desc_list[["GSE125856 : 1 vs 2"]]
      
      
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
          dplyr::select(
            -comparison, 
            -description,
            -Subtype_Immune_Model_Based,
            -cancer_type_abbreviation,
            -sample_type_id,
            -subtype_selected
          ) |>
          dplyr::filter(condition %in% local_condition) |> 
          dplyr::mutate(sample = gsub("-", "_", sample)) |> 
          tidyr::drop_na(age, gender, stage) |>
          dplyr::group_by(gender) |>
          dplyr::filter(dplyr::n() >= 5) |>
          dplyr::ungroup() |>
          dplyr::group_by(stage) |>
          dplyr::filter(dplyr::n() >= 5) |>
          dplyr::ungroup()
        
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
          covariates = covariate_arg,
          run_sva = TRUE
        )
      }, error = function(e) {
        message("paired_diff failed: ", 
                conditionMessage(e))
        NULL
      })
      
      
      ### Return data
      return(paireDiff)
      
    }  
  )
}

# Function to compute isoform fractions of diseases
## Takes as input a list of dataframes (one disease per dataframe with multiple studies)
disease_isoform_fractions <- function(dfs) {
  
  # For each dataframe do
  result <- do.call(rbind, lapply(seq_along(dfs), function(i) {
    
    df <- dfs[[i]]
    
    # Genes
    x <- df |>
      dplyr::select(.id, gene, padj_expression) |>
      dplyr::filter(padj_expression < 0.05) |>
      dplyr::distinct(.id, gene)
    
    # Isoforms
    y <- df |>
      dplyr::select(.id, gene, padj_splicing) |>
      dplyr::filter(padj_splicing < 0.05) |>
      dplyr::distinct(.id, gene)
    
    # n genes
    x_n <- x |>
      dplyr::count(.id, name = "n_expression")
    
    # n isoforms
    y_n <- y |>
      dplyr::count(.id, name = "n_splicing")
    
    # genes + isoforms
    union_n <- dplyr::bind_rows(
      x |> dplyr::select(.id, gene),
      y |> dplyr::select(.id, gene)
    ) |>
      dplyr::distinct(.id, gene) |>
      dplyr::count(.id, name = "n_union")
    
    # Isoform fraction
    fractions <- union_n |>
      dplyr::left_join(y_n, by = ".id") |>
      dplyr::mutate(
        n_splicing = tidyr::replace_na(n_splicing, 0),
        isoform_fraction = n_splicing / n_union
      )
    
    # Collecting all of the above in a dataframe
    data.frame(
      Disease = names(dfs)[i],
      Genes = median(x_n$n_expression, na.rm = TRUE),
      Isoforms = median(y_n$n_splicing, na.rm = TRUE),
      Genes_union_isoforms = median(union_n$n_union, na.rm = TRUE),
      Isoform_fraction = median(fractions$isoform_fraction, na.rm = TRUE)
    )
  }))
  
  result |>
    dplyr::arrange(Isoform_fraction)
}
# Function to run all tcga datasets with run_paired_diff
## Wrapper functions that loops over all tcga datasets curated with relevant run_paired_diff arguments
## At the end the object is saved as a .rds file format.
run_all_tcga <- function(
    tcga_dir, 
    n_cores, 
    out_dir) 
{
  # list all files 
  files <- list.files(tcga_dir, 
                      pattern = "\\.xlsx$", 
                      full.names = TRUE)
  
  # Create directory only if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  for (f in files) {
    file_name <- tools::file_path_sans_ext(basename(f))
    out_file <- file.path(out_dir, 
                          paste0(file_name, 
                                 "_paired_diff.rds"))
    
    # Skip if output file already exists
    if (file.exists(out_file)) {
      message("Skipping (already exists): ", 
              basename(out_file))
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
    
    readr::write_rds(paired_diff_results, 
                     file = out_file)
  }
}


# Function to run all archs4 datasets with run_paired_diff
## Wrapper functions that loops over all archs4 datasets curated with relevant run_paired_diff arguments.
## At the end the object is saved as a .rds file format.
run_all_archs4 <- function(
    datasets, 
    base_dir, 
    n_cores, 
    out_dir
) 
{
  # Create directory only if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  for (d in datasets) {
    
    sample_path <- file.path(base_dir, 
                             paste0(d, 
                                    "_samples.xlsx"))
    comparison_path <- file.path(base_dir, 
                                 paste0(d, 
                                        "_comparisons.xlsx"))
    
    out_file <- file.path(
      out_dir,
      paste0(d, "_paired_diff.rds")
    )
    
    # Skip if output already exists
    if (file.exists(out_file)) {
      message("Skipping (already exists): ", basename(out_file))
      next
    }
    
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
      file = file.path(out_dir, 
                       paste0(d, 
                              "_paired_diff.rds"))
    )
  }
}


# Function to combine cancer subtypes to one list of lists
combine_paired_diff_lists <- function(input_dir,
                                      output_dir,
                                      output_file = "combined_list.rds") {
  
  rds_files <- list.files(
    input_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )
  
  combined_list <- list()
  
  for (f in rds_files) {
    x <- readRDS(f)
    
    # append contents directly, not as combined_list[[file_name]]
    combined_list <- c(combined_list, x)
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(
    combined_list,
    file.path(output_dir, output_file)
  )
  
  invisible(combined_list)
}


# Function convert list from run_paired_diff to a tibble/dataframe
## and is used for any downstream analysis of run_paired_diff objects.
combine_lists_to_df <- function(
    pairedGSEA_list = paired_diff_results,
    archs4 = TRUE,
    n_cores
) 
{
  # Register n cores
  doMC::registerDoMC(cores = n_cores)
  if (archs4) {
    # apply tibble class to lists
    paired_tibble <- lapply(pairedGSEA_list, 
                            tibble::as_tibble)
    
    # convert each list to dataframe
    plyr::ldply(
      .data = paired_tibble,
      .parallel = TRUE
    )
  } else {
    paired_tibble <- as.data.frame(pairedGSEA_list@listData)
  }
  
  
}


# Function count expression and splicing features
## Functions to count DGE and DGU
counting_features <- function(
    path = "data/paired_diff_archs4/", 
    padj_score = "padj_splicing",
    new_path = "data/counting_features_archs4/",
    ...
)
{
  # Create new directory if not present
  if (!dir.exists(new_path)) {
    dir.create(new_path)
  }
  
  # List all files to perform counting of features on
  files <- list.files(path = path, 
                      pattern = "\\.rds$", 
                      full.names = TRUE)
  
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
                           paste0("feature_counts_", 
                                  padj_score, 
                                  ".rds"))
    
    # New path with addition of new filename
    out_path <- file.path(new_path, new_name)
    
    # Write new RDS file
    readr::write_rds(result, out_path)
    
  })
}


# Function count expression and splicing genesets
## Functions to count DGE and DGU (for genesets)
## Output is both each individually filtered padj score and the raw dataframe for each disease.
counting_genesets <- function(
    path_to_pd_objects,
    padj_score,
    out_dir,
    genesets = NULL,
    ...
) 
{
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
      combine_lists_to_df() |> 
      dplyr::mutate(gene = stringr::str_remove(gene, "\\..*"))
    
    # pairedGSEA over-representation analysis
    paired_ora_results <- plyr::dlply(
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
    
    # Build output file path
    out_file_raw_ora <- file.path(out_dir, 
                                  paste0(new_base_fname, 
                                         "ora_raw.rds"))
    
    # Save result
    saveRDS(paired_ora_results, 
            file = out_file_raw_ora)
    
    message("Saved: ", out_file_raw_ora)
    
    ora_results <- paired_ora_results |>
      dplyr::filter(.data[[padj_score]] < 0.05) |>
      dplyr::count(pathway, name = "n") |>
      dplyr::arrange(desc(n))
    
    # Build output file path
    out_file_count <- file.path(out_dir, paste0(new_base_fname, 
                                                "ora_results_", 
                                                padj_score ,
                                                ".rds"))
    
    # Save result
    saveRDS(ora_results, 
            file = out_file_count)
    
    message("Saved: ", 
            out_file_count)
  })
}

# ---- Preliminary RNA-seq Plots ----

# Function to extract disease base before generic filenames.
# Removes ".rds", "_paired_diff", and "_ora_raw"
extract_disease <- function(path) {
  basename(path) |>
    stringr::str_remove("\\.rds$") |>
    stringr::str_remove("_ora_raw$") |>
    stringr::str_remove("_paired_diff$")
}

# Summarize one file for either ARCHS4 and TCGA single ID diseases.
## Computes significant genes, isoforms, overlap, and only isoforms.
## ARCHS4 and TCGA is controlled by single_id. single_id = FALSE is ARCHS4 and single_id = TRUE is for TCGA
summarise_one_file <- function(
    path_to_file = path,
    single_id = FALSE,
    feature = gene,
    ...
) {
  feature <- rlang::enquo(feature)
  disease <- extract_disease(path_to_file)
  
  data_file <- readr::read_rds(path_to_file) |>
    combine_lists_to_df()
  
  # significant differential expression
  significant_genes <- data_file |>
    dplyr::filter(padj_expression < 0.05)
  
  # significant differential splicing
  significant_isoforms <- data_file |>
    dplyr::filter(padj_splicing < 0.05)
  
  if (!single_id) {
    
    # multi-.id mode
    diff_genes <- significant_genes |>
      dplyr::distinct(.id, 
                      !!feature) |>
      dplyr::count(.id, 
                   name = "Genes")
    
    diff_isoforms <- significant_isoforms |>
      dplyr::distinct(.id, 
                      !!feature) |>
      dplyr::count(.id, 
                   name = "Isoforms")
    
    diff_overlap <- dplyr::inner_join(
      significant_genes |> 
        dplyr::distinct(.id, 
                        !!feature),
      significant_isoforms |> 
        dplyr::distinct(.id, 
                        !!feature),
      by = c(".id", 
             "gene")
    ) |>
      dplyr::count(.id, 
                   name = "Overlap")
    
    diff_only_isoforms <- dplyr::anti_join(
      significant_isoforms |> 
        dplyr::distinct(.id, 
                        !!feature),
      significant_genes |> 
        dplyr::distinct(.id, 
                        !!feature),
      by = c(".id", 
             "gene")
    ) |>
      dplyr::count(.id, 
                   name = "Only_Isoforms")
    
    diff_list_features <- list(
      diff_isoforms,
      diff_genes,
      diff_overlap,
      diff_only_isoforms
    ) |>
      purrr::reduce(dplyr::full_join, 
                    by = ".id") |>
      tidyr::replace_na(list(
        Isoforms = 0,
        Genes = 0,
        Overlap = 0,
        Only_Isoforms = 0
      )) |>
      dplyr::mutate(disease = disease, 
                    .before = 1)
    
  } else {
    
    # single-.id mode
    # count everything in the file
    genes_tbl <- significant_genes |>
      dplyr::distinct(!!feature) |>
      dplyr::summarise(Genes = dplyr::n())
    
    isoforms_tbl <- significant_isoforms |>
      dplyr::distinct(!!feature) |>
      dplyr::summarise(Isoforms = dplyr::n())
    
    overlap_tbl <- dplyr::inner_join(
      significant_genes |> 
        dplyr::distinct(!!feature),
      significant_isoforms |> 
        dplyr::distinct(!!feature),
      by = !!feature
    ) |>
      dplyr::summarise(Overlap = dplyr::n())
    
    only_isoforms_tbl <- dplyr::anti_join(
      significant_isoforms |> 
        dplyr::distinct(!!feature),
      significant_genes |> 
        dplyr::distinct(!!feature),
      by = !!feature
    ) |>
      dplyr::summarise(Only_Isoforms = dplyr::n())
    
    diff_list_features <- dplyr::bind_cols(
      genes_tbl,
      isoforms_tbl,
      overlap_tbl,
      only_isoforms_tbl
    ) |>
      tidyr::replace_na(list(
        Isoforms = 0,
        Genes = 0,
        Overlap = 0,
        Only_Isoforms = 0
      )) |>
      dplyr::mutate(
        disease = disease,
        .id = disease,
        .before = 1
      )
  }
  
  diff_list_features
}

# summarize all files in a folder wrapper function. As the above function summarize one file at a this extend
## the function to all files. 
summarise_all_files <- function(input_dir,
                                pattern = "_paired_diff\\.rds$",
                                ...) {
  files <- list.files(
    path = input_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  purrr::map_dfr(files, summarise_one_file)
}

# Make disease-level medians per study per disease (used in boxplots)
medians_for_plot <- function(all_results,
                             group_col = "disease",
                             ...) {
  all_results |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::summarise(
      Isoforms = median(Isoforms, na.rm = TRUE),
      Genes = median(Genes, na.rm = TRUE),
      Overlap = median(Overlap, na.rm = TRUE),
      Only_Isoforms = median(Only_Isoforms, na.rm = TRUE),
      .groups = "drop"
    )
}

# Plotting significant features per study per disease 
plot_summary_feature_boxplots <- function(
    input_dir,
    pattern = "_paired_diff\\.rds$",
    output_dir = "plots",
    single_id = FALSE,
    width = 8,
    height = 7,
    dpi = 300,
    ...
) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir,
               recursive = TRUE, 
               showWarnings = FALSE)
  }
  
  # Running with diseases summarized
  all_results <- summarise_all_files(
    input_dir = input_dir,
    pattern = pattern,
    single_id = single_id
  )
  
  boxplot_df <- medians_for_plot(all_results, 
                                 group_col = "disease")
  
  disease_levels <- all_results |>
    dplyr::distinct(disease) |>
    dplyr::arrange(disease) |>
    dplyr::pull(disease)
  
  disease_colors <- viridisLite::viridis(
    n = length(disease_levels),
    option = "D"
  )
  
  names(disease_colors) <- disease_levels
  
  # use all_results for boxplots
  boxplot_long <- all_results |>
    tidyr::pivot_longer(
      cols = c(Isoforms, 
               Genes, 
               Overlap, 
               Only_Isoforms),
      names_to = "feature_type",
      values_to = "value"
    ) |>
    dplyr::mutate(
      value = dplyr::if_else(is.na(value),
                             0, 
                             value),
      value = value + 1
    )
  
  make_one_plot <- function(feature_name) {
    y_label <- dplyr::case_when(
      feature_name == "Isoforms" ~ "Significant Isoforms / Comparison",
      feature_name == "Genes" ~ "Significant Genes / Comparison",
      feature_name == "Overlap" ~ "Overlap Genes / Comparison",
      feature_name == "Only_Isoforms" ~ "Only Isoforms / Comparison",
      TRUE ~ feature_name
    )
    
    median_labels <- boxplot_long |>
      dplyr::filter(feature_type == feature_name) |>
      dplyr::group_by(disease) |>
      dplyr::summarise(
        median_value = median(value - 1, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        disease_label = paste0(
          disease,
          "\nMedian = ",
          median_value
        )
      )
    
    disease_axis_labels <- stats::setNames(
      median_labels$disease_label,
      median_labels$disease
    )
    
    p <- boxplot_long |>
      dplyr::filter(feature_type == feature_name) |>
      dplyr::mutate(
        disease = factor(disease,
                         levels = disease_levels),
        disease = forcats::fct_reorder(disease, 
                                       value, 
                                       .fun = median)
      ) |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = disease,
          y = value,
          color = disease,
          fill = disease
        )
      ) +
      ggbeeswarm::geom_beeswarm(
        dodge.width = 0.90,
        size = 2,
        alpha = 0.5
      ) +
      ggplot2::geom_boxplot(
        width = 0.7,
        color = "black",
        alpha = 0.5,
        outlier.shape = NA,
        lwd = 0.5
      ) +
      ggplot2::scale_x_discrete(
        labels = disease_axis_labels
      ) +
      ggplot2::scale_color_manual(
        values = disease_colors,
        drop = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = disease_colors,
        drop = FALSE
      ) +
      ggplot2::scale_y_log10(
        breaks = c(1, 
                   11, 
                   101, 
                   1001, 
                   10001),
        labels = c(0, 
                   10, 
                   100, 
                   1000, 
                   10000)
      ) +
      ggplot2::labs(
        x = "",
        y = y_label,
        title = paste0(
          "Significant features - ",
          stringr::str_replace_all(feature_name, 
                                   "_", 
                                   " ")
        )
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    ggplot2::ggsave(
      filename = file.path(output_dir, 
                           paste0("boxplot_", 
                                  feature_name, 
                                  ".png")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
    
    p
  }
  
  plots <- list(
    Isoforms = make_one_plot("Isoforms"),
    Genes = make_one_plot("Genes"),
    Overlap = make_one_plot("Overlap"),
    Only_Isoforms = make_one_plot("Only_Isoforms")
  )
  
  plots$All <- patchwork::wrap_plots(
    plots$Isoforms,
    plots$Genes,
    plots$Overlap,
    plots$Only_Isoforms,
    ncol = 2
  ) + 
    patchwork::plot_layout(
      guides = "collect",
      axis_titles = "collect"
    )
  
  ggplot2::ggsave(
    filename = file.path(output_dir, 
                         "boxplot_All.png"),
    plot = plots$All,
    width = width * 2,
    height = height * 2,
    dpi = dpi
  )
  
  return(list(
    all_results = all_results,
    boxplot_summary = boxplot_df,
    plots = plots
  ))
}


# ---- summarize one geneset dataframe .rds file ----
summarise_one_geneset_file <- function(path_to_file) {
  disease <- extract_disease(path_to_file)
  
  data_file <- readr::read_rds(path_to_file) |>
    tibble::as_tibble()
  
  significant_genes <- data_file |>
    dplyr::filter(padj_expression < 0.05)
  
  significant_isoforms <- data_file |>
    dplyr::filter(padj_splicing < 0.05)
  
  diff_genes <- significant_genes |>
    dplyr::group_by(.id) |>
    dplyr::count(.id, name = "Genesets")
  
  diff_isoforms <- significant_isoforms |>
    dplyr::group_by(.id) |>
    dplyr::count(.id, name = "Isoforms")
  
  diff_overlap <- dplyr::inner_join(
    significant_genes |> dplyr::distinct(.id, pathway),
    significant_isoforms |> dplyr::distinct(.id, pathway),
    by = c(".id", "pathway")
  ) |>
    dplyr::group_by(.id) |>
    dplyr::count(.id, name = "Overlap")
  
  overlap_paths <- dplyr::inner_join(
    significant_genes |> dplyr::distinct(.id, pathway),
    significant_isoforms |> dplyr::distinct(.id, pathway),
    by = c(".id", "pathway")
  )
  
  diff_only_isoforms <- dplyr::anti_join(
    significant_isoforms |> dplyr::select(.id, pathway),
    overlap_paths,
    by = c(".id", "pathway")
  ) |>
    dplyr::group_by(.id) |>
    dplyr::count(.id, name = "Only_Isoforms")
  
  list(
    diff_isoforms,
    diff_genes,
    diff_overlap,
    diff_only_isoforms
  ) |>
    purrr::reduce(dplyr::full_join, by = ".id") |>
    tidyr::replace_na(list(
      Isoforms = 0,
      Genesets = 0,
      Overlap = 0,
      Only_Isoforms = 0
    )) |>
    dplyr::mutate(disease = disease, .before = 1)
}


# ---- summarize all geneset files ----
summarise_all_geneset_files <- function(
    input_dir,
    pattern = "_ora_raw\\.rds$"
) {
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  purrr::map_dfr(files, summarise_one_geneset_file)
}


# ---- make disease-level medians ----
medians_for_plot <- function(all_results, 
                             group_col = "disease") {
  all_results |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::summarise(
      Isoforms = median(Isoforms, na.rm = TRUE),
      Genes = median(Genes, na.rm = TRUE),
      Overlap = median(Overlap, na.rm = TRUE),
      Only_Isoforms = median(Only_Isoforms, na.rm = TRUE),
      .groups = "drop"
    )
}


# ---- plot geneset boxplots ----
plot_summary_geneset_boxplots <- function(
    input_dir = "data/ora_archs4/",
    pattern = "_ora_raw\\.rds$",
    output_dir = "plots_genesets",
    width = 8,
    height = 7,
    dpi = 300
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_results <- summarise_all_geneset_files(
    input_dir = input_dir,
    pattern = pattern
  )
  
  boxplot_summary <- all_results |>
    dplyr::group_by(disease) |>
    dplyr::summarise(
      Isoforms = median(Isoforms, na.rm = TRUE),
      Genesets = median(Genesets, na.rm = TRUE),
      Overlap = median(Overlap, na.rm = TRUE),
      Only_Isoforms = median(Only_Isoforms, na.rm = TRUE),
      .groups = "drop"
    )
  
  boxplot_long <- all_results |>
    tidyr::pivot_longer(
      cols = c(Isoforms, Genesets, Overlap, Only_Isoforms),
      names_to = "feature_type",
      values_to = "value"
    ) |>
    dplyr::mutate(
      value_raw = dplyr::if_else(is.na(value), 0, value),
      value = value_raw + 1
    )
  
  make_one_plot <- function(feature_name) {
    y_label <- dplyr::case_when(
      feature_name == "Isoforms" ~ "Significant Splicing Genesets / Comparison",
      feature_name == "Genesets" ~ "Significant Expression Genesets / Comparison",
      feature_name == "Overlap" ~ "Overlap Genesets / Comparison",
      feature_name == "Only_Isoforms" ~ "Only Splicing Genesets / Comparison",
      TRUE ~ feature_name
    )
    
    disease_order <- boxplot_long |>
      dplyr::filter(feature_type == feature_name) |>
      dplyr::group_by(disease) |>
      dplyr::summarise(
        med = median(value_raw, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::arrange(med) |>
      dplyr::pull(disease)
    
    plot_data <- boxplot_long |>
      dplyr::filter(feature_type == feature_name) |>
      dplyr::mutate(
        disease = factor(disease, levels = disease_order)
      )
    
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = disease,
        y = value,
        color = disease,
        fill = disease
      )
    ) +
      ggplot2::geom_boxplot(
        width = 0.7,
        color = "black",
        alpha = 0.5,
        outlier.shape = NA,
        lwd = 0.5
      ) +
      ggbeeswarm::geom_beeswarm(
        dodge.width = 0.90,
        size = 2,
        alpha = 0.5
      ) +
      ggplot2::scale_y_log10(
        breaks = c(1, 11, 101, 1001, 10001),
        labels = c(0, 10, 100, 1000, 10000)
      ) +
      ggplot2::labs(
        x = "",
        y = y_label,
        title = paste0(
          "Significant genesets - ",
          stringr::str_replace_all(feature_name, "_", " ")
        )
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0("boxplot_", feature_name, ".png")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
    
    p
  }
  
  plots <- list(
    Isoforms = make_one_plot("Isoforms"),
    Genesets = make_one_plot("Genesets"),
    Overlap = make_one_plot("Overlap"),
    Only_Isoforms = make_one_plot("Only_Isoforms")
  )
  
  plots$All <- patchwork::wrap_plots(
    plots$Isoforms,
    plots$Genesets,
    plots$Overlap,
    plots$Only_Isoforms,
    ncol = 2
  ) +
    patchwork::plot_layout(
      guides = "collect",
      axis_titles = "collect"
    )
  
  ggplot2::ggsave(
    filename = file.path(output_dir, "boxplot_All.png"),
    plot = plots$All,
    width = width * 2,
    height = height * 2,
    dpi = dpi
  )
  
  list(
    all_results = all_results,
    boxplot_summary = boxplot_summary,
    boxplot_long = boxplot_long,
    plots = plots
  )
}



make_ridge_from_ora <- function(path, 
                                pattern = "_ora_raw\\.rds$", 
                                recursive = TRUE,
                                x_limits = c(-1, 1),
                                save_plot = TRUE,
                                plot_dir = "plots",
                                plot_name = "ridge_plot_faceted",
                                width = 16,
                                height = 7,
                                dpi = 300) {
  
  ora_files <- list.files(
    path = path,
    pattern = pattern,
    full.names = TRUE,
    recursive = recursive
  )
  
  if (length(ora_files) == 0) {
    stop("No files matching pattern found.")
  }
  
  filter_scenarios <- tibble::tribble(
    ~filter_type, ~filter_fun,
    "padj_expression < 0.05",
    \(df) dplyr::filter(df, padj_expression < 0.05),
    
    "padj_expression < 0.05 | padj_splicing < 0.05",
    \(df) dplyr::filter(df, padj_expression < 0.05 | padj_splicing < 0.05),
    
    "padj_splicing < 0.05",
    \(df) dplyr::filter(df, padj_splicing < 0.05)
  )
  
  all_ridge_data <- purrr::map_dfr(
    .x = stats::setNames(
      ora_files,
      stringr::str_remove(basename(ora_files), pattern)
    ),
    .f = function(file) {
      
      pairedOraResultsDF <- readRDS(file)
      
      purrr::map_dfr(
        seq_len(nrow(filter_scenarios)),
        function(i) {
          
          filtered_df <- filter_scenarios$filter_fun[[i]](pairedOraResultsDF) |>
            dplyr::group_by(.id) |>
            dplyr::filter(dplyr::n() > 10) |>
            dplyr::ungroup()
          
          filtered_df |> 
            dplyr::group_by(.id) |> 
            dplyr::summarize(
              corEnrichScores = stats::cor(
                enrichment_score_expression,
                enrichment_score_splicing,
                method = "spearman",
                use = "pairwise.complete.obs"
              ),
              n_sig_genesets = dplyr::n(),
              filter_type = filter_scenarios$filter_type[[i]],
              .groups = "drop"
            )
        }
      )
    },
    .id = "Disease"
  )
  
  disease_levels <- all_ridge_data |> 
    dplyr::group_by(Disease) |> 
    dplyr::summarize(
      overall_median = stats::median(corEnrichScores, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    dplyr::arrange(dplyr::desc(overall_median)) |> 
    dplyr::pull(Disease)
  
  all_ridge_data_plot <- all_ridge_data |> 
    dplyr::mutate(
      Disease = factor(Disease, levels = disease_levels),
      filter_type = factor(
        filter_type,
        levels = c(
          "padj_expression < 0.05",
          "padj_expression < 0.05 | padj_splicing < 0.05",
          "padj_splicing < 0.05"
        )
      )
    )
  
  medians_df_plot <- all_ridge_data_plot |> 
    dplyr::group_by(filter_type, Disease) |> 
    dplyr::summarize(
      median = stats::median(corEnrichScores, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    dplyr::mutate(
      Disease_num = as.numeric(Disease),
      median_label = paste0("median = ", round(median, 2))
    )
  
  ridge_plot <- all_ridge_data_plot |> 
    ggplot2::ggplot(ggplot2::aes(
      x = corEnrichScores,
      y = Disease,
      fill = Disease,
      color = Disease
    )) +
    ggridges::geom_density_ridges(
      alpha = 0.5,
      scale = 1.2,
      rel_min_height = 0.01
    ) +
    ggplot2::geom_segment(
      data = medians_df_plot,
      ggplot2::aes(
        x = median,
        xend = median,
        y = Disease_num,
        yend = Disease_num + 1,
        color = Disease
      ),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::geom_text(
      data = medians_df_plot,
      ggplot2::aes(
        x = median + 0.03,
        y = Disease_num + 0.5,
        label = median_label
      ),
      color = "black",
      inherit.aes = FALSE,
      hjust = 0,
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~ filter_type, nrow = 1) +
    ggplot2::coord_cartesian(xlim = x_limits, clip = "off") +
    ggplot2::labs(
      x = "Spearman's rho correlation coefficient",
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20),
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 14),
      plot.margin = ggplot2::margin(10, 40, 10, 10)
    )
  
  summary_table <- all_ridge_data |> 
    dplyr::group_by(filter_type, Disease) |> 
    dplyr::summarize(
      median = stats::median(corEnrichScores, na.rm = TRUE),
      n = sum(!is.na(corEnrichScores)),
      .groups = "drop"
    ) |> 
    dplyr::arrange(filter_type, median)
  
  if (save_plot) {
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
    
    file_path <- file.path(plot_dir, paste0(plot_name, ".png"))
    
    ggplot2::ggsave(
      filename = file_path,
      plot = ridge_plot,
      width = width,
      height = height,
      dpi = dpi
    )
    
    message("Plot saved to: ", file_path)
  }
  
  return(list(
    data = all_ridge_data,
    plot_data = all_ridge_data_plot,
    medians = medians_df_plot,
    plot = ridge_plot,
    summary = summary_table
  ))
}



# ---- Preliminary meta-analyses functions ----

background_cutoff <- function(
    dfs = list_of_dataframes,
    padj_thresh = 0.05,
    expression_padj_col = "padj_expression",
    splicing_padj_col = "padj_splicing",
    min_studies_range = 1:6,
    studies_comparison = c(">=", "=="),
    datasets_comparison = c(">=", "==")
) {
  
  studies_comparison  <- match.arg(studies_comparison)
  datasets_comparison <- match.arg(datasets_comparison)
  
  
  precompute_gene_counts <- function(padj_col) {
    dfs |>
      purrr::imap_dfr(~ .x |>
                        dplyr::filter(
                          !is.na(.data[[padj_col]]),
                          .data[[padj_col]] < padj_thresh
                        ) |>
                        dplyr::distinct(.id, gene) |>
                        dplyr::count(gene, name = "n_studies") |>
                        dplyr::mutate(dataset = .y)
      ) |>
      dplyr::group_by(gene) |>
      tidyr::nest(dataset_counts = c(dataset, n_studies)) |>
      dplyr::ungroup()
  }
  
  genes_for_threshold <- function(
    counts,
    min_studies,
    min_datasets,
    studies_comparison,
    datasets_comparison
  ) {
    
    counts |>
      dplyr::mutate(
        n_datasets = purrr::map_int(
          dataset_counts,
          ~ {
            if (studies_comparison == ">=") {
              sum(.x$n_studies >= min_studies)
            } else {
              sum(.x$n_studies == min_studies)
            }
          }
        )
      ) |>
      dplyr::filter(
        if (datasets_comparison == ">=") {
          n_datasets >= min_datasets
        } else {
          n_datasets == min_datasets
        }
      ) |>
      dplyr::pull(gene)
  }
  
  expression_counts <- precompute_gene_counts(expression_padj_col)
  splicing_counts   <- precompute_gene_counts(splicing_padj_col)
  
  threshold_results <- tidyr::expand_grid(
    min_studies = min_studies_range,
    min_datasets = seq_along(dfs)
  ) |>
    dplyr::mutate(
      expression_genes = purrr::map2(
        min_studies, min_datasets,
        ~ genes_for_threshold(
          expression_counts,
          .x,
          .y,
          studies_comparison,
          datasets_comparison
        )
      ),
      splicing_genes = purrr::map2(
        min_studies, min_datasets,
        ~ genes_for_threshold(
          splicing_counts,
          .x,
          .y,
          studies_comparison,
          datasets_comparison
        )
      ),
      n_expression_genes = purrr::map_int(expression_genes, length),
      n_splicing_genes = purrr::map_int(splicing_genes, length),
      n_overlap_genes = purrr::map2_int(
        expression_genes,
        splicing_genes,
        ~ length(intersect(.x, .y))
      )
    ) |>
    dplyr::select(-expression_genes, -splicing_genes)
  
  threshold_results_long <- threshold_results |>
    tidyr::pivot_longer(
      cols = c(n_expression_genes, n_splicing_genes, n_overlap_genes),
      names_to = "signal_type",
      values_to = "n_genes"
    ) |>
    dplyr::mutate(
      signal_type = dplyr::recode(
        signal_type,
        n_expression_genes = "Expression",
        n_splicing_genes = "Splicing",
        n_overlap_genes = "Overlap"
      ),
      signal_type = factor(
        signal_type,
        levels = c("Expression", "Splicing", "Overlap")
      )
    )
  
  plot <- ggplot2::ggplot(
    threshold_results_long,
    ggplot2::aes(
      x = min_studies,
      y = min_datasets,
      fill = n_genes
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = n_genes),
      color = "white",
      size = 5,
      angle = -45
    ) +
    ggplot2::facet_wrap(~ signal_type) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(threshold_results_long$min_studies))
    ) +
    ggplot2::scale_y_continuous(
      breaks = sort(unique(threshold_results_long$min_datasets))
    ) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 12),
      strip.text = ggplot2::element_text(size = 14)
    ) +
    ggplot2::labs(
      x = "Minimum studies",
      y = "Minimum diseases",
      fill = "Genes retained"
    )
  
  list(
    threshold_results = threshold_results,
    threshold_results_long = threshold_results_long,
    plot = plot
  )
}


sampled_cutoff <- function(
    dfs,
    padj_thresh = 0.05,
    n_repeats = 99,
    seed = 123,
    expression_padj_col = "padj_expression",
    splicing_padj_col = "padj_splicing",
    feature_col = "gene",
    min_studies_range = 1:6
) {
  set.seed(seed)
  
  test_thresholds <- tidyr::expand_grid(
    min_studies = min_studies_range,
    min_datasets = seq_along(dfs)
  )
  
  full_universe <- dfs |>
    purrr::imap_dfr(~ .x |>
                      dplyr::distinct(
                        dataset = .y,
                        .id,
                        feature = .data[[feature_col]]
                      )
    )
  
  feature_universe <- full_universe |>
    dplyr::distinct(dataset, feature)
  
  get_x_by_dataset <- function(padj_col) {
    dfs |>
      purrr::imap_dfr(~ .x |>
                        dplyr::filter(
                          !is.na(.data[[padj_col]]),
                          .data[[padj_col]] < padj_thresh
                        ) |>
                        dplyr::distinct(dataset = .y, feature = .data[[feature_col]])
      ) |>
      dplyr::count(dataset, name = "n_to_sample")
  }
  
  run_signal <- function(padj_col, signal_type) {
    x_by_dataset <- get_x_by_dataset(padj_col)
    
    purrr::map_dfr(seq_len(n_repeats), function(repeat_id) {
      
      sampled_features <- feature_universe |>
        dplyr::left_join(x_by_dataset, by = "dataset") |>
        dplyr::mutate(n_to_sample = tidyr::replace_na(n_to_sample, 0L)) |>
        dplyr::group_split(dataset) |>
        purrr::map_dfr(function(df) {
          n_to_sample <- unique(df$n_to_sample)
          n_to_sample <- min(n_to_sample, nrow(df))
          
          df |>
            dplyr::slice_sample(n = n_to_sample)
        }) |>
        dplyr::select(dataset, feature)
      
      sampled_summary <- full_universe |>
        dplyr::semi_join(sampled_features, by = c("dataset", "feature")) |>
        dplyr::count(dataset, feature, name = "n_studies") |>
        dplyr::count(feature, n_studies, name = "n_datasets")
      
      test_thresholds |>
        dplyr::mutate(
          n_genes = purrr::map2_int(
            min_studies,
            min_datasets,
            ~ sampled_summary |>
              dplyr::filter(
                n_studies == .x,
                n_datasets == .y
              ) |>
              nrow()
          ),
          repeat_id = repeat_id,
          signal_type = signal_type
        )
    })
  }
  
  repeat_results <- dplyr::bind_rows(
    run_signal(expression_padj_col, "expression"),
    run_signal(splicing_padj_col, "splicing")
  )
  
  threshold_results_median <- repeat_results |>
    dplyr::group_by(signal_type, min_studies, min_datasets) |>
    dplyr::summarise(
      median_n_genes = stats::median(n_genes),
      mean_n_genes = mean(n_genes),
      sd_n_genes = stats::sd(n_genes),
      .groups = "drop"
    )
  
  plot <- ggplot2::ggplot(
    threshold_results_median,
    ggplot2::aes(
      x = min_studies,
      y = min_datasets,
      fill = median_n_genes
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = median_n_genes), 
                       color = "white", 
                       size = 4,
                       angle = -45
                       ) +
    ggplot2::facet_wrap(~ signal_type) +
    ggplot2::scale_x_continuous(breaks = sort(unique(threshold_results_median$min_studies))) +
    ggplot2::scale_y_continuous(breaks = sort(unique(threshold_results_median$min_datasets))) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Minimum studies",
      y = "Minimum diseases",
      fill = "Median random features retained"
    )
  
  list(
    repeat_results = repeat_results,
    threshold_results_median = threshold_results_median,
    plot = plot
  )
}


enrich_cutoff <- function(
    observed_results,
    sampled_summary
) {
  
  observed_long <- observed_results |>
    dplyr::select(
      min_studies,
      min_datasets,
      n_expression_genes,
      n_splicing_genes,
      n_overlap_genes
    ) |>
    tidyr::pivot_longer(
      cols = c(
        n_expression_genes,
        n_splicing_genes,
        n_overlap_genes
      ),
      names_to = "signal_type",
      values_to = "observed"
    ) |>
    dplyr::mutate(
      signal_type = dplyr::recode(
        signal_type,
        n_expression_genes = "expression",
        n_splicing_genes = "splicing",
        n_overlap_genes = "overlap"
      )
    )
  
  enrichment_results <- sampled_summary |>
    dplyr::rename(
      expected_mean = mean_n_genes,
      expected_median = median_n_genes,
      expected_sd = sd_n_genes
    ) |>
    dplyr::left_join(
      observed_long,
      by = c(
        "signal_type",
        "min_studies",
        "min_datasets"
      )
    ) |>
    dplyr::mutate(
      
      enrichment =
        observed / expected_median,
      
      log2_enrichment =
        log2((observed) / (expected_median)),
      
      z_score =
        dplyr::if_else(
          expected_sd > 0,
          (observed - expected_median) / expected_sd,
          NA_real_
        ),
      
      empirical_p =
        NA_real_
    )
  
  enrichment_results
}

extract_background_genes <- function(
    dfs,
    padj_col,
    min_studies,
    min_datasets,
    padj_thresh = 0.05,
    feature_col = "gene",
    cumulative = TRUE,
    new_filename = NULL
) {
  
  feature_counts <- dfs |>
    purrr::imap_dfr(~ .x |>
                      dplyr::filter(
                        !is.na(.data[[padj_col]]),
                        .data[[padj_col]] < padj_thresh
                      ) |>
                      dplyr::distinct(
                        dataset = .y,
                        .id,
                        feature = .data[[feature_col]]
                      )
    ) |>
    dplyr::count(dataset, feature, name = "n_studies")
  
  result <- if (cumulative) {
    
    feature_counts |>
      dplyr::filter(n_studies >= min_studies) |>
      dplyr::count(feature, name = "n_datasets") |>
      dplyr::filter(n_datasets >= min_datasets) |>
      dplyr::pull(feature)
    
  } else {
    
    feature_counts |>
      dplyr::filter(n_studies == min_studies) |>
      dplyr::count(feature, name = "n_datasets") |>
      dplyr::filter(n_datasets == min_datasets) |>
      dplyr::pull(feature)
  }
  
  if (!is.null(new_filename)) {
    saveRDS(result, file = new_filename)
  }
  
  result
}

get_disease_patterns <- function(
    dfs,
    padj_col,
    min_studies = 6,
    min_datasets = 9,
    padj_thresh = 0.05,
    feature_col = "gene",
    cumulative = TRUE
) {
  
  feature_counts <- dfs |>
    purrr::imap_dfr(~ .x |>
                      dplyr::filter(
                        !is.na(.data[[padj_col]]),
                        .data[[padj_col]] < padj_thresh
                      ) |>
                      dplyr::distinct(
                        dataset = .y,
                        .id,
                        feature = .data[[feature_col]]
                      )
    ) |>
    dplyr::count(dataset, feature, name = "n_studies")
  
  if (cumulative) {
    feature_counts |>
      dplyr::group_by(feature) |>
      dplyr::filter(sum(n_studies >= min_studies) >= min_datasets) |>
      dplyr::arrange(dataset, .by_group = TRUE) |>
      dplyr::summarise(
        diseases = paste(dataset, collapse = " | "),
        study_counts = paste(n_studies, collapse = " | "),
        combo = paste0(dataset, ":", n_studies, collapse = " | "),
        n_datasets_passing = sum(n_studies >= min_studies),
        .groups = "drop"
      ) |>
      dplyr::rename(!!feature_col := feature)
  } else {
    feature_counts |>
      dplyr::group_by(feature) |>
      dplyr::filter(sum(n_studies == min_studies) == min_datasets) |>
      dplyr::arrange(dataset, .by_group = TRUE) |>
      dplyr::summarise(
        diseases = paste(dataset, collapse = " | "),
        study_counts = paste(n_studies, collapse = " | "),
        combo = paste0(dataset, ":", n_studies, collapse = " | "),
        n_datasets_passing = sum(n_studies == min_studies),
        .groups = "drop"
      ) |>
      dplyr::rename(!!feature_col := feature)
  }
}

# Run fgsea::fora
run_fgsea_on_cutoff <- function(dfs = dfs,
                                path_to_feature_cutoff) {
  
  pathways <- readr::read_rds("data/genesets/genesets.rds")
  features <- readr::read_rds(path_to_feature_cutoff)
  universe_features <- dplyr::bind_rows(dfs) |>
    dplyr::filter(if_any(starts_with("padj"), ~ !is.na(.))) |>
    dplyr::distinct(gene, .keep_all = TRUE) |>
    dplyr::pull(gene)
  
  fgsea::fora(pathways,
              genes = features,
              universe = universe_features,
              minSize = 25) |> 
    dplyr::filter(padj < 0.05)
}

# Prep function for geneSetSimplifyR
combine_fgsea_gene_isoform <- function(gene_res = fgsea_cumu_d9_s6_genes,
                                       isoform_res = fgsea_cumu_d9_s6_isoforms,
                                       sigGenes_genes = readr::read_rds("data/cut_off_cumulative/genes_diseases9_studies6.rds"),
                                       sigGenes_isoforms = readr::read_rds("data/cut_off_cumulative/genes_diseases9_studies6.rds"),
                                       save_list_path,
                                       save_enrichment_path,
                                       archs4_file = "/home/databases/archs4/v2.latest/human_tpm_v2.latest.h5",
                                       padj_cutoff = 0.05,
                                       pseudocount = 0.06) {
  
  universe_genes <- rhdf5::h5read(
    file = archs4_file,
    name = "meta/transcripts/ensembl_gene"
  ) |>
    unique()
  
  universe_transcripts <- rhdf5::h5read(
    file = archs4_file,
    name = "meta/transcripts/ensembl_id"
  )
  
  compute_enrichment <- function(res,
                                 sigFeatures,
                                 universe_features,
                                 feature_type,
                                 source) {
    
    n_sig <- length(sigFeatures)
    n_universe <- length(universe_features)
    
    tibble::as_tibble(res) |>
      dplyr::filter(padj < padj_cutoff) |>
      dplyr::mutate(
        feature_type = feature_type,
        source = source,
        enrichment_score = log2(
          ((overlap / size) / (n_sig / n_universe)) + pseudocount
        )
      )
  }
  
  combined_df <- dplyr::bind_rows(
    compute_enrichment(
      res = gene_res,
      sigFeatures = sigGenes_genes,
      universe_features = universe_genes,
      feature_type = "genes",
      source = "Expression"
    ),
    compute_enrichment(
      res = isoform_res,
      sigFeatures = sigGenes_isoforms,
      universe_features = universe_transcripts,
      feature_type = "isoforms",
      source = "Splicing"
    )
  )
  
  combined_lists <- combined_df |>
    dplyr::filter(!is.na(pathway)) |>
    dplyr::select(pathway, overlapGenes) |>
    dplyr::mutate(
      overlapGenes = purrr::map(overlapGenes, \(x) unique(as.character(x)))
    ) |>
    tibble::deframe()
  
  saveRDS(combined_lists, file = save_list_path)
  
  exampleEnrichment <- combined_df |>
    dplyr::select(source, pathway, enrichment_score)
  
  saveRDS(exampleEnrichment, file = save_enrichment_path)
  
  return(combined_df)
}


plot_volcano_expression <- function(
    df_list,
    gene_col = "gene",
    lfc_col = "lfc_expression",
    pval_col = "pvalue_expression",
    pval_threshold = 0.05,
    lfc_threshold = 1,
    label_top_n = 10
) {
  
  volcano_plots <- purrr::imap(df_list, function(df, name) {
    
    df_plot <- df |>
      dplyr::group_by(.data[[gene_col]]) |>
      dplyr::summarise(
        lfc_expression = stats::median(.data[[lfc_col]], na.rm = TRUE),
        pvalue_expression = stats::median(.data[[pval_col]], na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        neg_log10_pval = -log10(pvalue_expression),
        regulation = dplyr::case_when(
          pvalue_expression < pval_threshold &
            lfc_expression > lfc_threshold ~ "Up-regulated gene",
          pvalue_expression < pval_threshold &
            lfc_expression < -lfc_threshold ~ "Down-regulated gene",
          TRUE ~ "Not significant"
        )
      )
    
    top_labels <- df_plot |>
      dplyr::filter(regulation %in% c("Up-regulated gene", "Down-regulated gene")) |>
      dplyr::group_by(regulation) |>
      dplyr::slice_max(
        order_by = abs(lfc_expression),
        n = label_top_n,
        with_ties = FALSE
      ) |>
      dplyr::ungroup()
    
    ggplot2::ggplot(
      df_plot,
      ggplot2::aes(
        x = lfc_expression,
        y = neg_log10_pval,
        color = regulation
      )
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 2) +
      ggrepel::geom_text_repel(
        data = top_labels,
        ggplot2::aes(label = .data[[gene_col]]),
        size = 3,
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      ggplot2::geom_vline(
        xintercept = c(-lfc_threshold, lfc_threshold),
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(pval_threshold),
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::scale_color_manual(values = c(
        "Up-regulated gene" = "#F8766D",
        "Down-regulated gene" = "#619CFF",
        "Not significant" = "gray"
      )) +
      ggplot2::labs(
        title = paste("Expression Volcano Plot:", name),
        x = "lfc_expression",
        y = "-log10 pvalue_expression",
        color = "Expression regulation"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        text = ggplot2::element_text(size = 12)
      )
  })
  
  return(volcano_plots)
}

plot_splicing_scatter <- function(
    df_list,
    gene_col = "gene",
    dif_col = "max_abs_dif_splicing",
    pval_col = "pvalue_splicing",
    pval_threshold = 0.05,
    dif_threshold = 0.1,
    n_labels = 20
) {
  
  scatter_plots <- purrr::imap(df_list, function(df, name) {
    
    df_plot <- df |>
      dplyr::group_by(.data[[gene_col]]) |>
      dplyr::summarise(
        dif = stats::median(.data[[dif_col]], na.rm = TRUE),
        pval = stats::median(.data[[pval_col]], na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        neg_log10_pval = -log10(pval),
        regulation = dplyr::case_when(
          pval < pval_threshold &
            dif > dif_threshold ~ "Significant splicing change",
          TRUE ~ "Not significant"
        )
      )
    
    labels_df <- df_plot |>
      dplyr::filter(regulation != "Not significant") |>
      dplyr::slice_min(
        order_by = pval,
        n = n_labels,
        with_ties = FALSE
      )
    
    ggplot2::ggplot(
      df_plot,
      ggplot2::aes(
        x = dif,
        y = neg_log10_pval,
        color = regulation
      )
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 2) +
      ggrepel::geom_text_repel(
        data = labels_df,
        ggplot2::aes(label = .data[[gene_col]]),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.5,
        point.padding = 0.3,
        show.legend = FALSE
      ) +
      ggplot2::geom_vline(
        xintercept = dif_threshold,
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(pval_threshold),
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::scale_color_manual(values = c(
        "Significant splicing change" = "#F8766D",
        "Not significant" = "gray"
      )) +
      ggplot2::labs(
        title = paste("Splicing Scatter Plot:", name),
        x = dif_col,
        y = paste0("-log10(", pval_col, ")"),
        color = "Splicing significance"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        text = ggplot2::element_text(size = 12)
      )
  })
  
  return(scatter_plots)
}
# ---- Formal p-value integration functions ----

meta_pvalues_by_gene <- function(disease_dfs, min_p = 2e-16) {
  
  fix_p <- function(p) {
    p[p == 0] <- min_p
    p
  }
  
  combine_sumlog <- function(p) {
    p <- p |>
      stats::na.omit() |>
      as.numeric() |>
      fix_p()
    
    if (length(p) == 0) return(NA_real_)
    if (length(p) == 1) return(p)
    
    metap::sumlog(p)$p
  }
  
  compute_se <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }
  
  disease_dfs |>
    purrr::map(\(df) {
      df |>
        dplyr::mutate(
          dplyr::across(
            c(pvalue_expression, pvalue_splicing),
            fix_p
          )
        ) |>
        dplyr::summarise(
          expression_meta_p = combine_sumlog(pvalue_expression),
          splicing_meta_p   = combine_sumlog(pvalue_splicing),
          
          n_expression = sum(!is.na(pvalue_expression)),
          n_splicing   = sum(!is.na(pvalue_splicing)),
          
          median_lfc_expression = median(lfc_expression, na.rm = TRUE),
          median_lfc_splicing   = median(lfc_splicing, na.rm = TRUE),
          
          median_max_abs_dif_splicing =
            median(max_abs_dif_splicing, na.rm = TRUE),
          
          se_weighted_expression = compute_se(lfc_expression),
          se_weighted_splicing   = compute_se(lfc_splicing),
          
          .by = gene
        ) |>
        dplyr::mutate(
          expression_meta_fdr = stats::p.adjust(expression_meta_p, method = "BH"),
          splicing_meta_fdr   = stats::p.adjust(splicing_meta_p, method = "BH")
        )
    })
}


plot_volcano_meta <- function(
    df_list,
    gene_col = "gene",
    x_metric = "weighted_lfc_expression",
    padj_col = "expression_meta_fdr",
    padj_threshold = 0.05,
    log2fc_threshold = 1,
    min_p = 2e-16
) {
  
  volcano_plots <- purrr::imap(df_list, function(df, name) {
    
    df_plot <- df |>
      dplyr::mutate(
        neg_log10_padj = -log10(pmax(.data[[padj_col]], min_p)),
        regulation = dplyr::case_when(
          .data[[padj_col]] < padj_threshold &
            .data[[x_metric]] > log2fc_threshold ~ "Up-regulated",
          
          .data[[padj_col]] < padj_threshold &
            .data[[x_metric]] < -log2fc_threshold ~ "Down-regulated",
          
          TRUE ~ "Not significant"
        )
      )
    
    labels_df <- dplyr::bind_rows(
      df_plot |>
        dplyr::filter(regulation == "Up-regulated") |>
        dplyr::slice_min(order_by = .data[[padj_col]], n = 20, with_ties = FALSE),
      
      df_plot |>
        dplyr::filter(regulation == "Down-regulated") |>
        dplyr::slice_min(order_by = .data[[padj_col]], n = 20, with_ties = FALSE)
    )
    
    ggplot2::ggplot(
      df_plot,
      ggplot2::aes(
        x = .data[[x_metric]],
        y = neg_log10_padj,
        color = regulation
      )
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 2) +
      ggrepel::geom_text_repel(
        data = labels_df,
        ggplot2::aes(label = .data[[gene_col]]),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.2,
        show.legend = FALSE
      ) +
      ggplot2::geom_vline(
        xintercept = c(-log2fc_threshold, log2fc_threshold),
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(padj_threshold),
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::scale_color_manual(values = c(
        "Up-regulated" = "#F8766D",
        "Down-regulated" = "#619CFF",
        "Not significant" = "gray"
      )) +
      ggplot2::labs(
        title = paste("Volcano Plot:", name),
        x = x_metric,
        y = "-log10 adjusted p-value",
        color = "Regulation"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        text = ggplot2::element_text(size = 12)
      )
  })
  
  return(volcano_plots)
}

meta_pvalues_across_diseases <- function(meta_dfs, min_p = 2e-16) {
  
  fix_p <- function(p) {
    p[p == 0] <- min_p
    p
  }
  
  combine_p <- function(p, method = c("fisher", "edgington")) {
    method <- match.arg(method)
    
    p <- p |>
      stats::na.omit() |>
      as.numeric() |>
      fix_p()
    
    if (length(p) == 0) return(NA_real_)
    if (length(p) == 1) return(p)
    
    switch(
      method,
      fisher    = metap::sumlog(p)$p,
      edgington = metap::sump(p)$p
    )
  }
  
  weighted_lfc <- function(lfc, se) {
    keep <- !is.na(lfc) & !is.na(se) & se > 0
    
    if (!any(keep)) {
      return(median(lfc, na.rm = TRUE))
    }
    
    w <- 1 / se[keep]^2
    
    sum(w * lfc[keep]) / sum(w)
  }
  
  weighted_lfc_se <- function(se) {
    se <- se[!is.na(se) & se > 0]
    
    if (length(se) == 0) return(NA_real_)
    
    w <- 1 / se^2
    
    sqrt(1 / sum(w))
  }
  
  meta_dfs |>
    dplyr::bind_rows(.id = "disease") |>
    dplyr::summarise(
      expression_fisher_p =
        combine_p(expression_meta_p, "fisher"),
      expression_edgington_p =
        combine_p(expression_meta_p, "edgington"),
      
      splicing_fisher_p =
        combine_p(splicing_meta_p, "fisher"),
      splicing_edgington_p =
        combine_p(splicing_meta_p, "edgington"),
      
      n_diseases_expression =
        sum(!is.na(expression_meta_p)),
      n_diseases_splicing =
        sum(!is.na(splicing_meta_p)),
      
      weighted_lfc_expression =
        weighted_lfc(
          weighted_lfc_expression,
          se_weighted_expression
        ),
      
      weighted_lfc_splicing =
        weighted_lfc(
          weighted_lfc_splicing,
          se_weighted_splicing
        ),
      
      se_weighted_expression =
        weighted_lfc_se(se_weighted_expression),
      
      se_weighted_splicing =
        weighted_lfc_se(se_weighted_splicing),
      
      .by = gene
    ) |>
    dplyr::mutate(
      expression_fisher_fdr =
        stats::p.adjust(expression_fisher_p, method = "BH"),
      
      expression_edgington_fdr =
        stats::p.adjust(expression_edgington_p, method = "BH"),
      
      splicing_fisher_fdr =
        stats::p.adjust(splicing_fisher_p, method = "BH"),
      
      splicing_edgington_fdr =
        stats::p.adjust(splicing_edgington_p, method = "BH")
    )
}


run_fgsea_multilevel <- function(counted_ranks,
                                 pathways,
                                 scoreType = "pos",
                                 eps = 0) {
  
  fgsea_results <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = counted_ranks,
    scoreType = scoreType,
    eps = eps
  )
  
  return(data.frame(fgsea_results))
}


# ----  Cross-Indication Target Discovery ----
run_chembl_targets <- function(data,
                               method = c("fora", "fgseaMultilevel"),
                               chembl,
                               gene_map,
                               padj_cutoff = 0.05) {
  
  method <- match.arg(method)
  
  #----------------------------------------------------------
  # Helper: map ENSG IDs to HGNC symbols
  #----------------------------------------------------------
  map_ensg_to_symbol <- function(genes, gene_map) {
    tibble::tibble(gene_raw = genes) |>
      dplyr::mutate(
        ensembl_gene_id = stringr::str_remove(gene_raw, "\\.\\d+$")
      ) |>
      dplyr::left_join(
        gene_map,
        by = "ensembl_gene_id"
      ) |>
      dplyr::mutate(
        gene_symbol = dplyr::if_else(
          is.na(hgnc_symbol) | hgnc_symbol == "",
          ensembl_gene_id,
          hgnc_symbol
        )
      ) |>
      dplyr::pull(gene_symbol)
  }
  
  #----------------------------------------------------------
  # Prepare targets
  #----------------------------------------------------------
  if (method == "fora") {
    
    data <- data |>
      dplyr::mutate(
        overlapGenes = purrr::map(
          overlapGenes,
          map_ensg_to_symbol,
          gene_map = gene_map
        )
      )
    
    targets <- data |>
      dplyr::filter(padj < padj_cutoff) |>
      dplyr::select(pathway, padj, overlapGenes) |>
      tidyr::unnest_longer(overlapGenes) |>
      dplyr::rename(gene = overlapGenes) |>
      dplyr::mutate(
        source = "fora",
        NES = NA_real_,
        direction = "unknown"
      )
    
  } else {
    
    data <- data |>
      dplyr::mutate(
        leadingEdge = purrr::map(
          leadingEdge,
          map_ensg_to_symbol,
          gene_map = gene_map
        )
      )
    
    targets <- data |>
      dplyr::filter(padj < padj_cutoff) |>
      dplyr::select(pathway, padj, NES, leadingEdge) |>
      tidyr::unnest_longer(leadingEdge) |>
      dplyr::rename(gene = leadingEdge) |>
      dplyr::mutate(
        source = "fgseaMultilevel",
        direction = dplyr::case_when(
          NES > 0 ~ "up_enriched",
          NES < 0 ~ "down_enriched",
          TRUE ~ "neutral"
        )
      )
  }
  
  #----------------------------------------------------------
  # ChEMBL join
  #----------------------------------------------------------
  chembl <- chembl |>
    dplyr::rename(gene = hgnc_symbol)
  
  hits <- targets |>
    dplyr::filter(!is.na(gene), gene != "") |>
    dplyr::inner_join(
      chembl,
      by = "gene"
    )
  
  #----------------------------------------------------------
  # Compound ranking
  #----------------------------------------------------------
  compound_scores <- hits |>
    dplyr::group_by(compound_chembl_id) |>
    dplyr::summarise(
      max_phase = dplyr::first(phase),
      n_target_genes = dplyr::n_distinct(gene),
      n_pathways = dplyr::n_distinct(pathway),
      target_genes = paste(sort(unique(gene)), collapse = ";"),
      pathways = paste(sort(unique(pathway)), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(n_pathways),
      dplyr::desc(n_target_genes)
    )
  
  list(
    targets = targets,
    hits = hits,
    compound_scores = compound_scores
  )
}
