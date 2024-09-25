source("config.r")


if (!("pheatmap" %in% installed.packages())) {
  install.packages("pheatmap")
}
library(pheatmap)

#Heatmap generation for the top 50 genes

top_fifty_genes <- top_fifty$Gene
top_fifty_expression_df <- filtered_expression_df[top_fifty_genes, ]

sample_annotations <- data.frame(
  Group = filtered_metadata$time_status
)

rownames(sample_annotations) <- colnames(scaled_expression_df)

scaled_expression_df <- t(scale(t(top_fifty_expression_df)))

heatmap_plot <- pheatmap(
  scaled_expression_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = sample_annotations,
  main = "Heatmap of significant genes"
)

ggsave(filename = file.path(plots_dir, "SRP062829P_heatmop_top_50_expressed_genes.png"),
       plot = heatmap_plot,)

#Heatmap for all genes

top_genes <- deseq_df$Gene
top_expression_df <- filtered_expression_df[top_genes, ]

significant_genes <- deseq_df %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  dplyr::arrange(padj)

significant_gene_names <- significant_genes$Gene

significant_expression_df <- filtered_expression_df[significant_gene_names, ]

scaled_significant_expression_df <- t(scale(t(significant_expression_df)))

sample_annotations <- data.frame(
  Group = filtered_metadata$time_status
)

rownames(sample_annotations) <- colnames(scaled_expression_df)

scaled_expression_df <- t(scale(t(top_expression_df)))

heatmap_plot_top_50 <- pheatmap(
  significant_expression_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = sample_annotations,
  main = "Heatmap of Top 50 Significant Genes"
)

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

ggsave(filename = file.path(plots_dir, "SRP062829P_heatmop_top_50_expressed_genes.png"),
       plot = heatmap_plot_top_50,)
