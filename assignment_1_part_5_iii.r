source("config.r")

if(!("gprofiler2" %in% installed.packages())) {
  BiocManager::install("gprofiler2", update = FALSE)
}

# Load required libraries
library(DESeq2)
library(dplyr)
library(gprofiler2)
library(readr)
library(tibble)

# Load your data files (adjust file paths if needed)
filtered_data_file <- "results/filtered_data.tsv"
filtered_metadata_file <- "results/filtered_metadata.tsv"

# Read the filtered data (gene expression counts)
filtered_expression_df <- readr::read_tsv(filtered_data_file) %>%
  tibble::column_to_rownames("Symbol")

# Read the filtered metadata (experimental design information)
filtered_metadata <- readr::read_tsv(filtered_metadata_file)

# Convert 'time_status' in metadata to a factor (this variable will be used for design)
filtered_metadata <- filtered_metadata %>%
  dplyr::mutate(time_status = as.factor(time_status))

# Create DESeq2 object from the filtered data
ddset <- DESeqDataSetFromMatrix(
  countData = filtered_expression_df,  # The gene expression counts data
  colData = filtered_metadata,         # The metadata (experimental design info)
  design = ~ time_status               # Design formula using 'time_status'
)

# Run DESeq2 analysis to get differential expression results
deseq_object <- DESeq(ddset)

# Extract the results from the DESeq2 analysis
deseq_results <- results(deseq_object)

# Shrink the log2 fold changes (optional, but improves results)
deseq_results <- lfcShrink(
  deseq_object, 
  coef = 2,  # Coefficient for log2 fold change shrinkage (depends on your design)
  res = deseq_results
)

# Convert DESeq2 results to a data frame for easier handling
deseq_df <- deseq_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::mutate(threshold = padj < 0.05) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange))

# Save the full differential expression results as a TSV file (optional)
readr::write_tsv(
  deseq_df,
  file.path(results_dir, "gProfiler2_Gene_Ontology_differential_expression_results.tsv")
)

# Extract the list of significantly differentially expressed genes (adjusted p-value < 0.05)
significant_genes <- deseq_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::pull(Gene)

# Run Gene Set Enrichment Analysis (GSEA) using gProfiler2 and Gene Ontology (GO)
gsea_results <- gost(
  query = significant_genes,    # List of significant genes
  organism = "hsapiens",         # Specify the organism (use "hsapiens" for human)
  sources = "GO",                # Use "GO" for Gene Ontology (you can change this to other ontologies)
  user_threshold = 0.05          # Significance threshold for enrichment
)

# Convert the GSEA results to a data frame
gsea_df <- as.data.frame(gsea_results$result)

# Save the full GSEA results to a TSV file
readr::write_tsv(
  gsea_df,
  file.path(results_dir, "gProfiler2_Gene_Ontology_gsea_go_enrichment_results.tsv")
)

# Extract the top 20 enriched terms by p-value for easier reporting
top_enrichments <- gsea_df %>%
  dplyr::arrange(p_value) %>%
  dplyr::slice(1:20)

# Save the top 20 enriched terms to a TSV file
readr::write_tsv(
  top_enrichments,
  file.path(results_dir, "gProfiler2_Gene_Ontology_top_20_gsea_go_enrichment_results.tsv")
)

# Summary: You now have full differential expression results and GSEA results.
