source("config.r")

# make sure the differential expression results file exists
dge_mapped_df_file <- file.path(results_dir, "SRP062829_diff_expr_results.tsv")
if (!(file.exists(dge_mapped_df_file))) {
    source("assignment_1_part_3.r")
}

# install packages

if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}

if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}

# load libraries
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(magrittr)

# set seed
set.seed(123)

# read in data

dge_mapped_df <- readr::read_tsv(
  file.path(results_dir, "SRP062829_diff_expr_results.tsv")
) %>%
  dplyr::select(Gene, log2FoldChange, padj)

# create a list of gene sets
mm_hallmark_sets <- msigdbr(
  species = "Mus musculus", 
  category = "H"
)

# Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- dge_mapped_df$log2FoldChange
names(lfc_vector) <- dge_mapped_df$Gene

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

# Convert the GSEA results to a data frame
gsea_result_df <- data.frame(gsea_results@result)

# Write the GSEA results to a TSV file
readr::write_tsv(
  gsea_result_df,
  file.path(
    results_dir,
    "SRP062829_clusterProfiler_HallmarkONTO_gsea_results.tsv"
  )
)

gsea_top_ten <- gsea_result_df %>%
  dplyr::slice_max(NES, n = 10)

# Write the top 10 GSEA results to a TSV file
readr::write_tsv(
  gsea_top_ten,
  file.path(
    results_dir,
    "SRP062829_clusterProfiler_HallmarkONTO_top_10_gsea_results.tsv"
  )
)