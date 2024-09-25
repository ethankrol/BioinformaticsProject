source("config.r")

mapped_expression_df_file <- file.path(results_dir, "mapped_expression_df.tsv")

# Install the Mouse package
if (!("org.Mm.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}

if (!("tidyr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("tidyr")
}

# Load the required libraries
library(org.Mm.eg.db)
library(dplyr)
library(magrittr)
library(tidyr)

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

mapped_list <- mapIds(
  org.Mm.eg.db, # Replace with annotation package for your organism
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "list"
)
# Map Ensembl IDs to their associated ??? can't be hugo cuz their not human
# TO DO: find out what gene ids to map to
mapped_df <- mapped_list %>%
    tibble::enframe(name = "ENSEMBL", value = "SYMBOL") %>%
      # enframe() makes a `list` column; we will simplify it with unnest()
      # This will result in one row of our data frame per list item
      tidyr::unnest(cols = SYMBOL)

multi_mapped <- mapped_df %>%
  # Let's count the number of times each Ensembl ID appears in `Ensembl` column
  dplyr::count(ENSEMBL, name = "symbol_id_count") %>%
  # Arrange by the genes with the highest number of Entrez IDs mapped
  dplyr::arrange(desc(symbol_id_count))

collapsed_mapped_df <- mapped_df %>%
  # Group by Ensembl IDs
  dplyr::group_by(ENSEMBL) %>%
  # Collapse the Entrez IDs `mapped_df` into one column named `all_entrez_ids`
  dplyr::summarize(all_symbol_ids = paste(SYMBOL, collapse = ";"))

collapsed_mapped_df %>%
  # Filter `collapsed_mapped_df` to include only the rows where
  # `all_entrez_ids` values include the ";" character --
  # these are the rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_symbol_ids, ";")) %>%
  # We only need a preview here
  head()

final_mapped_df <- data.frame(
  "first_mapped_symbol_id" = mapIds(
    org.Mm.eg.db, # Replace with annotation package for your organism
    keys = expression_df$Gene,
    keytype = "ENSEMBL", # Replace with the gene identifiers used in your data
    column = "SYMBOL", # The type of gene identifiers you would like to map to
    multiVals = "first" # Keep only the first mapped value for each Ensembl ID
  )
) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("ENSEMBL") %>%
  # Add the multiple mappings data from `collapsed_mapped_df` using Ensembl IDs
  dplyr::inner_join(collapsed_mapped_df, by = "ENSEMBL") %>%
  # Now let's add on the rest of the expression data
  dplyr::inner_join(expression_df, by = c("ENSEMBL" = "Gene"))

final_mapped_df %>%
  # Filter `final_mapped_df` to rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_symbol_ids, ";")) %>%
  head()

readr::write_tsv(final_mapped_df, file.path(
  results_dir,
  "SRP062829_Symbol_IDs.tsv" # Replace with a relevant output file name
))

mapped_data <- readr::read_tsv("results/SRP062829_Symbol_IDs.tsv")
