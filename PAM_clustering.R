# Load necessary libraries
library(cluster)
library(dplyr)
library(ggplot2)
library(readr)

results_dir <- "results_assignment3"
num_clusters <- 5

# Function to perform PAM clustering from scratch (relying on pam())
perform_pam_clustering <- function(data, num_clusters) {
  # Scale the data (mean = 0, variance = 1)
  data_scaled <- scale(data)
  
  # Run PAM clustering
  pam_result <- pam(data_scaled, k = num_clusters)
  
  # Extract results manually
  medoids <- pam_result$medoids          # Medoids
  clustering <- pam_result$clustering    # Cluster assignments for each point
  sil_widths <- pam_result$silinfo$widths # Silhouette widths
  
  # Display results
  cat("Medoids:\n")
  print(medoids)
  cat("\nCluster Assignments:\n")
  print(clustering)

  # Visualization: PCA Plot
  pca_res <- prcomp(data_scaled)
  pca_df <- as.data.frame(pca_res$x)
  pca_df$Cluster <- as.factor(clustering)

  ggplot(pca_df, aes(PC1, PC2, color = Cluster)) +
    geom_point(size = 3) +
    labs(title = paste0("PCA Plot for ", num_clusters, " clusters")) +
    theme_minimal()

  # Visualization: Silhouette Plot
  plot(silhouette(clustering, dist(data_scaled)),
       main = paste0("Silhouette Plot for ", num_clusters, " clusters"))
  
  return(pam_result)
}

num_genes_list <- c(10, 100, 1000, 5000, 10000)

for (num_genes in num_genes_list) {
  expression_data <- read_tsv(file.path(results_dir, paste0("top_", num_genes, "_variable_genes.tsv")), col_names = TRUE)
  
  
  # Separate the 'Symbol' column (gene identifiers) from the numeric data
  expression_numeric <- expression_data %>% select(-Symbol)  # Keep only numeric columns
  
  # Transpose the data so that each row is a sample
  expression_numeric <- t(expression_numeric)
  
  # Convert transposed matrix to a data frame and assign sample names as row names
  expression_numeric <- as.data.frame(expression_numeric)
  rownames(expression_numeric) <- colnames(expression_data)[-1]  # Set row names to sample names
  
  # Scale the data (mean = 0, variance = 1)
  expression_scaled <- scale(expression_numeric)
  
  # Perform PAM clustering for chosen amount of clusters
  pam_result <- perform_pam_clustering(expression_scaled, num_clusters)
  
  clusters <- pam_result$clustering
  
  # Load metadata
  metadata <- read_tsv("results/filtered_metadata.tsv", col_names = TRUE)
  
  # Add cluster info to metadata
  clusters_df <- data.frame(refinebio_accession_code = rownames(expression_numeric), Cluster = clusters)
  metadata_with_clusters <- left_join(metadata, clusters_df, by = "refinebio_accession_code")
  
  # Replace NAs with empty strings
  metadata_with_clusters[is.na(metadata_with_clusters)] <- ""
  
  # Write out the results
  write_tsv(metadata_with_clusters, file.path(results_dir, paste0(num_clusters, "_pam_clusters_from_", num_genes, "_genes.tsv")))
}

