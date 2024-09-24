source("config.r")

## install packages
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}

## load libraries
library(DESeq2)
library(ggplot2)
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
    stringr::str_equal(refinebio_time, "24h post accelerated defeat") ~ "Shortly_After_Defeat",
    stringr::str_equal(refinebio_time, "48h post csds") ~ "Shortly_After_Defeat",
    stringr::str_equal(refinebio_time, "28d post csds") ~ "One_Month_After_Defeat",
    stringr::str_equal(refinebio_time, "28d + 1h stress post csds") ~ "One_Month_Post_Defeat_+_1_Hour_Stress",
    TRUE ~ "NA"
  )) %>%
  dplyr::mutate(time_status = factor(time_status))

filtered_metadata <- metadata %>%
  dplyr::filter(time_status != "NA") %>%
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

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = filtered_metadata,
  # Supply our experimental variable to `design`
  design = ~time_status
)

# Run the DESeq2 analysis
deseq_object <- DESeq(ddset)

deseq_results <- results(deseq_object)

deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

## convert to dataframe
# this is of class DESeqResults -- we want a data frame
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))


# write to tsv
readr::write_tsv(
  deseq_df,
  file.path(
    results_dir,
    "SRP062829_diff_expr_results.tsv" 
  )
)

# Create a volcano plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "SRP062829_volcano_plot.png")
)