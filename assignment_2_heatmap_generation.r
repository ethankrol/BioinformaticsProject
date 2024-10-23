library(mclust)
library(pheatmap)

filtered_expression_df <- readr::read_tsv(/results/filtered_data.tsv) %>%
  tibble::column_to_rownames("Symbol")
gene_variances <- apply(filtered_expression_df, 1, var)
top_5000_genes <- order(gene_variances, decreasing = TRUE)[1:5000]
top_5000_gene_data <- filtered_expression_df[top_5000_genes, ]

#Generate MClust clustering for 5000 genes if needed
mclust_top_5000 = t(top_5000_gene_data)
data_scaled <- scale(mclust_top_5000)
gmm_model_5000 = Mclust(data_scaled)
clusters_5000 <- gmm_model_5000$classification

scaled_top_5000_genes <- t(scale(t(top_5000_gene_data)))

filtered_metadata <- readr::read_tsv("/results/filtered_metadata.tsv")

#These are the row annotations for the heatmap. The variables are just a column with the clustering groups for each of the 55 samples.

heatmap_annotations <- data.frame(
  Cluster_5000_genes = as.factor(clusters_5000),
  Sample_Group= as.factor(filtered_metadata$time_status)
)


rownames(heatmap_annotations) <- colnames(data_5000)

heatmap_plot <- pheatmap(
  scaled_top_5000_genes,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_col = heatmap_annotations,
  main = "Heatmap of 5,000 Most Variable Genes with Clustering Annotations"
)
