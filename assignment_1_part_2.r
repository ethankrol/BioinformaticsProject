source("config.r")

if (!exists("filtered_data_file") | !exists("filtered_metadata_file")) {
  source("generate_filtered_data.r")
}

if (!file.exists(filtered_data_file)| !file.exists(filtered_metadata_file)) {
  source("generate_filtered_data.r")
}

## load libraries
library(ggplot2)

filtered_expression_df <- readr::read_tsv(filtered_data_file) %>%
  tibble::column_to_rownames("Symbol")
filtered_metadata <- readr::read_tsv(filtered_metadata_file)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = filtered_expression_df,
  # Supply the `colData` with our metadata data frame
  colData = filtered_metadata,
  # Supply our experimental variable to `design`
  design = ~time_status
)
head(filtered_expression_df)

colnames(colData(ddset))
vsd <- vst(ddset, blind = FALSE)
plotPCA(vsd, intgroup = c("time_status", "refinebio_subject"))
