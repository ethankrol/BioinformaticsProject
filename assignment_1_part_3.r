source("config.r")

if (!exists("filtered_data_file") | !exists("filtered_metadata_file")) {
  source("generate_filtered_data.r")
}

if (!file.exists(filtered_data_file)| !file.exists(filtered_metadata_file)) {
  source("generate_filtered_data.r")
}

## install packages
if (!("EnhancedVolcano" %in% installed.packages())) {
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("apeglm", update = FALSE)
}

## load libraries
library(ggplot2)
library(dplyr)
library(DESeq2)

set.seed(123)

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

# create a table of top fifty differentially expressed genes
top_fifty <- deseq_df %>%
    dplyr::slice(1:50)

# write to tsv
readr::write_tsv(
  top_fifty,
  file.path(
    results_dir,
    "SRP062829_top_fifty_diffexpr_genes.tsv" 
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