source("config.r")

## install packages
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("pheatmap" %in% installed.packages())) {
  BiocManager::install("pheatmap", update = FALSE)
}

## load libraries
library(DESeq2)
library(pheatmap)
library(magrittr)
library(ggplot2)

set.seed(123)

# read in metadata
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

# round all expression counts
filtered_expression_df <- round(filtered_expression_df)

# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = filtered_expression_df, # the counts values for all samples
  colData = metadata, # annotation data for the samples
  design = ~1 # Here we are not specifying a model
  # Replace with an appropriate design variable for your analysis
)