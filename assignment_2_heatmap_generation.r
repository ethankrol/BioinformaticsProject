library(pheatmap)

scaled_top_5000_genes <- t(scale(t(top_5000_gene_data)))

filtered_metadata <- readr::read_tsv("/results/filtered_metadata.tsv")

spectral_data <- readr::read_tsv("/results_assignment3/2_spectral_clusters_from_5000_genes.tsv")

Kmeans_data <- readr::read_tsv("/results_assignment3/K_means/K_means_6_clusters_from_5000_genes.tsv")

pam_data <- readr::read_tsv("/results_assignment3/5_pam_clusters_from_5000_genes.tsv")

affinity_data <- readr::read_tsv("results_assignment3/9_affinity_propagation_clusters_from_5000_genes.tsv")

heatmap_annotations <- data.frame(
  GMM_Cluster_5000_genes = as.factor(clusters_5000),
  K_means_cluster_5000_genes = as.factor(Kmeans_data$Cluster),
  spectral_cluster_5000_genes = as.factor(spectral_data$Cluster),
  pam_cluster_5000_genes = as.factor(pam_data$Cluster),
  affinity_cluster_5000_genes = as.factor(affinity_data$Cluster),
  Sample_Group= as.factor(filtered_metadata$time_status)
)

ann_colors <- list(
  GMM_Cluster_5000_genes = c("1" = "blue", "2" = "cyan", "3" = "magenta", "4" = "orange", "5" = "green", "6" = "purple", "7" = "pink"),
  K_means_cluster_5000_genes = c("0" = "red", "1" = "yellow", "2" = "limegreen", "3" = "skyblue", "4" = "salmon", "5" = "gold"),
  spectral_cluster_5000_genes = c("0" = "darkblue", "1" = "lightgreen"),
  pam_cluster_5000_genes = c("1" = "darkred", "2" = "darkgreen", "3" = "darkorange", "4" = "darkviolet", "5" = "gray"),
  affinity_cluster_5000_genes = c("0" = "lightpink", "1" = "lightcyan", "2" = "orange", "3" = "lightgreen", "4" = "lightblue", "5" = "lightcoral", "6" = "lightgoldenrod", "7" = "lightsalmon", "8" = "lightseagreen"),
  Sample_Group = c("One_Day_After_Defeat" = "red", "One_Month_After_Defeat" = "blue")  
)


rownames(heatmap_annotations) <- colnames(data_5000)

heatmap_plot <- pheatmap(
  scaled_top_5000_genes,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_col = heatmap_annotations,
  annotation_colors = ann_colors,
  main = "Heatmap of 5,000 Most Variable Genes with Clustering Annotations"
)
