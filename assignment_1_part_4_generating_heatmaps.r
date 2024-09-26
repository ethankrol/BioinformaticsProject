source("config.r")


if (!("pheatmap" %in% installed.packages())) {
  install.packages("pheatmap")
}
library(dplyr)
library(readr)
library(tibble)

#Heatmap generation for the top 50 genes


#Read in necessary files
top_fifty <- readr::read_tsv("results/SRP062829_top_fifty_diffexpr_genes.tsv")

filtered_expression_df <- readr::read_tsv(filtered_data_file) %>%
  tibble::column_to_rownames("Symbol")

filtered_metadata <- readr::read_tsv(filtered_metadata_file)

#Take the top fifty genes and normalize them
top_fifty_genes <- top_fifty$Gene
top_fifty_expression_df <- filtered_expression_df[top_fifty_genes, ]

scaled_expression_top_50_df <- t(scale(t(top_fifty_expression_df)))

#Create a dataframe to render the side bar
sample_annotations <- data.frame(
  Group = filtered_metadata$time_status
)

rownames(sample_annotations) <- colnames(scaled_expression_top_50_df)

#Plot the heatmap

heatmap_plot <- pheatmap(
  scaled_expression_top_50_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = sample_annotations,
  main = "Heatmap of top 50 significant genes"
)

#Save the plot
ggsave(filename = file.path(plots_dir, "SRP062829P_heatmop_top_50_expressed_genes.png"),
       plot = heatmap_plot_top_50,)

#Heatmap for all genes

#Identify genes with an appropriate p value and log-fold change
significant_genes <- deseq_df %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  dplyr::arrange(padj)

#Do the same process for the top 50 genes but with all of the significant genes
significant_gene_names <- significant_genes$Gene

significant_expression_df <- filtered_expression_df[significant_gene_names, ]

scaled_significant_expression_df <- t(scale(t(significant_expression_df)))

heatmap_plot <- pheatmap(
  scaled_significant_expression_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = sample_annotations,
  main = "Heatmap of significant genes"
)

ggsave(filename = file.path(plots_dir, "SRP062829P_heatmop_top_expressed_genes.png"),
       plot = heatmap_plot,)
