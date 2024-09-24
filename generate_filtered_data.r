## install packages
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}

## load libraries
library(DESeq2)
library(magrittr)

set.seed(123)

metadata <- readr::read_tsv(metadata_file)
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

metadata <- metadata %>%
  # Let's get the RPL10 mutation status from this variable
  dplyr::mutate(time_status = dplyr::case_when(
    stringr::str_equal(refinebio_time, "24h post accelerated defeat") ~ "One_Day_After_Defeat",
    stringr::str_equal(refinebio_time, "48h post csds") ~ "Two_Days_After_Defeat",
    stringr::str_equal(refinebio_time, "28d post csds") ~ "One_Month_After_Defeat",
    stringr::str_equal(refinebio_time, "28d + 1h stress post csds") ~ "One_Month_Post_Defeat_+_1_Hour_Stress",
    TRUE ~ "NA"
  )) %>%
  dplyr::mutate(time_status = factor(time_status))

filtered_metadata <- metadata %>%
  dplyr::filter(time_status != "NA") %>%
  dplyr::filter(time_status != "Two_Days_After_Defeat") %>%
  dplyr::filter(time_status != "One_Month_Post_Defeat_+_1_Hour_Stress")

# Filter the expression data to only include the samples in the metadata
filtered_expression_df <- expression_df %>%
  dplyr::select(all_of(filtered_metadata$refinebio_accession_code))

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- filtered_expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

# round all expression counts
gene_matrix <- round(filtered_expression_df)

# write the filtered data to a file
filtered_data_file <- file.path(results_dir, "filtered_data.tsv")
filtered_metadata_file <- file.path(results_dir, "filtered_metadata.tsv")
readr::write_tsv(gene_matrix, file.path(filtered_data_file))
readr::write_tsv(filtered_metadata, file.path(filtered_metadata_file))