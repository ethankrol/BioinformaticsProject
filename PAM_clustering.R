# Load necessary libraries
library(cluster)
library(dplyr)
library(ggplot2)
library(readr)
library(networkD3)


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

# Initialize vectors to store chi-squared test results
p_values <- c()
adjusted_p_values <- c()
test_names <- c()

# Loop through each pair of gene sets
for (i in 1:(length(num_genes_list) - 1)) {
  for (j in (i + 1):length(num_genes_list)) {
    
    num_genes1 <- num_genes_list[i]
    num_genes2 <- num_genes_list[j]
    
    # Read the clustering results for the two gene sets
    df1 <- read_tsv(paste0("results_assignment3/K_means/K_means_6_clusters_from_", num_genes1, "_genes.tsv"))
    df2 <- read_tsv(paste0("results_assignment3/K_means/K_means_6_clusters_from_", num_genes2, "_genes.tsv"))
    
    # Create a contingency table comparing the cluster assignments for the two sets
    contingency_table <- table(df1$Cluster, df2$Cluster)
    
    # Perform chi-squared test
    chi_test <- chisq.test(contingency_table)
    
    # Store p-value and test name
    p_values <- c(p_values, chi_test$p.value)
    test_names <- c(test_names, paste0(num_genes1, " genes vs ", num_genes2, " genes"))
  }
}

# Adjust p-values for multiple hypothesis testing using the Benjamini-Hochberg method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Combine results into a data frame
chi_squared_results <- data.frame(
  Test = test_names,
  Original_p_value = p_values,
  Adjusted_p_value = adjusted_p_values
)

# Print the chi-squared test results
print(chi_squared_results)

# Optionally, write the results to a file
write_tsv(chi_squared_results, "results_assignment3/chi_squared_results.tsv")

source <- c()
target <- c()
value <- c()

# Loop through each pair of gene sets
for (i in 1:(length(num_genes_list) - 1)) {
  num_genes1 <- num_genes_list[i]
  num_genes2 <- num_genes_list[i + 1]
  
  # Read the data for the two gene sets
  df1 <- read_tsv(paste0("results_assignment3/5_PAM_clusters_from_", num_genes1, "_genes.tsv"))
  df2 <- read_tsv(paste0("results_assignment3/5_PAM_clusters_from_", num_genes2, "_genes.tsv"))
  
  # Create a contingency table (cross-tabulation) between clusters from two different gene sets
  contingency_table <- table(df1$Cluster, df2$Cluster)
  
  # Extract the number of connections between clusters
  for (j in 1:nrow(contingency_table)) {
    for (k in 1:ncol(contingency_table)) {
      # Source cluster (in the first gene set) is row index j
      # Target cluster (in the second gene set) is column index k
      source <- c(source, 2 * (i - 1) + (j - 1))
      target <- c(target, 2 * (i - 1) + 2 + (k - 1))
      value <- c(value, as.integer(contingency_table[j, k]))
    }
  }
}

# Define node labels
node_labels <- unlist(lapply(num_genes_list, function(num_genes) {
  paste0(num_genes, " genes<br>Cluster ", 1:6)
}))

# Create Sankey diagram data
sankey_data <- list(
  nodes = data.frame(name = node_labels),
  links = data.frame(source = source, target = target, value = value)
)

# Create the Sankey diagram using networkD3
sankeyNetwork(Links = sankey_data$links, Nodes = sankey_data$nodes,
              Source = "source", Target = "target", Value = "value", NodeID = "name",
              fontSize = 12, nodeWidth = 30)