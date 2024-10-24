library(mclust)
library(ggplot2)
library(ggalluvial)
#Since our genes are currently our rows and our samples are currently our columns, we need to transpose this for mclust

mclust_top_5000 = t(top_5000_gene_data)
mclust_top_10000 = t(top_10000_gene_data)
mclust_top_1000 = t(top_1000_gene_data)
mclust_top_100 = t(top_100_gene_data)
mclust_top_10 = t(top_10_gene_data)

data_scaled <- scale(mclust_top_5000)
data_10000_scaled <- scale(mclust_top_10000)
data_1000_scaled <- scale(mclust_top_1000)
data_100_scaled <- scale(mclust_top_100)
data_10_scaled <- scale(mclust_top_10)

head(data_scaled)

#Create MClust model based off normalized data
gmm_model_5000 = Mclust(data_scaled)

#Find number of clusters
num_clusters <- gmm_model_5000$G
print(num_clusters)

#Find the classification of each sample in each cluster
clusters_5000 <- gmm_model_5000$classification
print(clusters_5000)
plot(gmm_model_5000, what = 'classification', dimens = c(40,29))

gmm_model_10000 = Mclust(data_10000_scaled)
num_clusters <- gmm_model_10000$G
print(num_clusters)
clusters_10000 <- gmm_model_10000$classification
print(clusters_10000)
plot(gmm_model_10000, what = 'classification', dimens = c(40,29))

gmm_model_1000 = Mclust(data_1000_scaled)
num_clusters <- gmm_model_1000$G
print(num_clusters)
clusters_1000 <- gmm_model_1000$classification
print(clusters_1000)
plot(gmm_model_1000, what = 'classification', dimens = c(40,29))

gmm_model_100 = Mclust(data_100_scaled)
num_clusters <- gmm_model_100$G
print(num_clusters)
clusters_100 <- gmm_model_100$classification
print(clusters_100)
plot(gmm_model_100, what = 'classification', dimens = c(40,29))

gmm_model_10 = Mclust(data_10_scaled)
num_clusters <- gmm_model_10$G
print(num_clusters)
clusters_10 <- gmm_model_10$classification
print(clusters_10)
plot(gmm_model_10, what = 'classification', dimens = c(1,2))

#Generate all combinations of gene cluster results
clustering_results <- list(
  "10 Genes" = clusters_10,
  "100 Genes" = clusters_100,
  "1000 Genes" = clusters_1000,
  "5000 Genes" = clusters_5000,
  "10000 Genes" = clusters_10000
)

chi_sq_results <- list()

combinations <- combn(names(clustering_results), 2, simplify = FALSE)

#For each pair, compute chi squared test on clustering results, and adjust p-value.
for (pair in combinations) {
  contingency_table <- table(clustering_results[[pair[1]]], clustering_results[[pair[2]]])
  
  chi_sq_test <- chisq.test(contingency_table)
  
  chi_sq_results[[paste(pair[1], "vs", pair[2], sep = " ")]] <- chi_sq_test$p.value
}

chi_sq_results_df <- data.frame(
  Comparison = names(chi_sq_results),
  P_value = unlist(chi_sq_results)
)

write.csv(chi_sq_results_df, "chi_squared_results_original_clusters.csv", row.names = FALSE)

chi_sq_results_df$Adjusted_P_value <- p.adjust(chi_sq_results_df$P_value)

#Transform the data sets so that each variable is a column. We need each column
# on the graph to be the number of genes fed into each clustering algorithm, and
#the y axis to be the different samples.

clusters_10_df = as.data.frame(clusters_10)

clusters_df <- data.frame(
  Sample = 1:length(clusters_10),
  Clusters_10_Genes = as.factor(clusters_10),
  Clusters_100_Genes = as.factor(clusters_100),
  Clusters_1000_Genes = as.factor(clusters_1000),
  Clusters_5000_Genes = as.factor(clusters_5000),
  Clusters_10000_Genes = as.factor(clusters_10000)
)

metadata_10_genes <- filtered_metadata
metadata_10_genes$cluster <- clusters_10

metadata_100_genes <- filtered_metadata
metadata_100_genes$cluster <- clusters_100

metadata_1000_genes <- filtered_metadata
metadata_1000_genes$cluster <- clusters_1000

metadata_5000_genes <- filtered_metadata
metadata_5000_genes$cluster <- clusters_5000

metadata_10000_genes <- filtered_metadata
metadata_10000_genes$cluster <- clusters_10000

write.csv(metadata_10_genes, "GMM_10_genes.csv")
write.csv(metadata_100_genes, "GMM_100_genes.csv")
write.csv(metadata_1000_genes, "GMM_1000_genes.csv")
write.csv(metadata_5000_genes, "GMM_5000_genes.csv")
write.csv(metadata_10000_genes, "GMM_10000_genes.csv")

ggplot(clusters_df,
       aes(y = 1, axis1 = Clusters_10_Genes, axis2 = Clusters_100_Genes, axis3 = Clusters_1000_Genes, axis4 = Clusters_5000_Genes, axis5 = Clusters_10000_Genes)) +
  geom_alluvium(aes(fill = Sample)) +
  geom_stratum(width = 1/6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes"),
                   expand = c(0.05, 0.05)) +
  labs(title = "Sankey Plot of Clusters Differing Across Numbers of Genes",
       x = "Number of Genes", y = "Samples") +
  theme_minimal()


