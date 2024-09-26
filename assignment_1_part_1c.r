source("config.r")

if (!exists("filtered_data_file") | !exists("filtered_metadata_file")) {
  source("generate_filtered_data.r")
}

if (!file.exists(filtered_data_file)| !file.exists(filtered_metadata_file)) {
  source("generate_filtered_data.r")
}

library(readr)
library(dplyr)
library(ggplot2)

# read in metadata
metadata <- readr::read_tsv(filtered_metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(filtered_data_file) %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Symbol")

# display dimensions
dim(expression_df)

# display first few rows
head(expression_df)

# count the number of genes
num_genes <- ncol(expression_df)
cat("Number of samples:", num_genes, "\n")

# count the number of samples
num_samples <- nrow(expression_df)
cat("Number of genes:", num_samples, "\n")

# log-scale the data
log_expression_df <- log2(expression_df + 1)

# calculate the median expression of each gene
gene_median_expression <- apply(log_expression_df, 1, median)
# Convert the vector into a data frame
gene_medians_df <- data.frame(Median = gene_median_expression)

# density plot of median expression
ggplot(gene_medians_df, aes(x = Median)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density plot of median gene expression",
       x = "Median gene expression",
       y = "Density")

# save the plot
ggsave(file.path(plots_dir, "median_expression_density_plot.png"))