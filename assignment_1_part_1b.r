source("config.r")

# Install the Mouse package
if (!("org.Mm.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}

# Load the required libraries
library(org.Mm.eg.db)
library(magrittr)

# read in metadata
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Bring back the "Gene" column in preparation for mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

# Map Ensembl IDs to their associated ??? can't be hugo cuz their not human
# TO DO: find out what gene ids to map to
map_df <- data.frame( 
  mapIds(
    org.Mm.eg.db, # Replace with annotation package for your organism
    keys = expression_df$Gene,
    keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
    column = "SYMBOL", # The type of gene identifiers you would like to map to
    multiVals = "first"
  )
) %>%
  dplyr::rename("Symbol" = "mapIds.org.Mm.eg.db..keys...expression_df.Gene..keytype....ENSEMBL...") %>%
  tibble::rownames_to_column("ENSEMBL")

mapped_expression_df <- expression_df %>%
  inner_join(map_df, by = c("Gene" = "ENSEMBL")) %>%
  dplyr::select(Symbol, everything()) %>%
  dplyr::select(-Gene)

# write to tsv
readr::write_tsv(mapped_expression_df, file.path(results_dir, "mapped_expression_df.tsv"))
