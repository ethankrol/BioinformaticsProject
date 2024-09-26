
source("config.r")

if (!exists("filtered_data_file") | !exists("filtered_metadata_file")) {
  source("generate_filtered_data.r")
}

if (!file.exists(filtered_data_file)| !file.exists(filtered_metadata_file)) {
  source("generate_filtered_data.r")
}

if (!("Rtsne" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("Rtsne", update = FALSE)
}

if (!("umap" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("umap", update = FALSE)
}

## load libraries
library(ggplot2)
library(dplyr)
library(DESeq2)

filtered_expression_df <- readr::read_tsv(filtered_data_file) %>%
  tibble::column_to_rownames("Symbol")
filtered_metadata <- readr::read_tsv(filtered_metadata_file)

# Convert 'time_status' to a factor to avoid warning
filtered_metadata <- filtered_metadata %>%
  dplyr::mutate(time_status = as.factor(time_status))

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = filtered_expression_df,
  # Supply the `colData` with our metadata data frame
  colData = filtered_metadata,
  # Supply our experimental variable to `design`
  design = ~time_status
)
head(filtered_expression_df)

colnames(colData(ddset))
vsd <- vst(ddset, blind = FALSE)
plotPCA(vsd, intgroup = c("time_status", "refinebio_subject"))

# Save the plot
ggsave(file.path(plots_dir, "pca_plot_1.png"))

library(Rtsne)

# Extract the variance-stabilized data matrix for t-SNE
vsd_matrix <- assay(vsd)

tsne_out <- Rtsne(t(vsd_matrix), dims = 2, perplexity = 15)

# Create a data frame for ggplot
tsne_df <- data.frame(
  X = tsne_out$Y[, 1],
  Y = tsne_out$Y[, 2],
  time_status = filtered_metadata$time_status
)

# Plot t-SNE results
ggplot(tsne_df, aes(x = X, y = Y, color = time_status)) +
  geom_point(size = 2) +
  ggtitle("t-SNE Plot of Samples")

# Save the plot
ggsave(file.path(plots_dir, "tsne_plot.png"))

library(umap)

# Run UMAP on the transposed data (samples as rows)
umap_out <- umap::umap(t(vsd_matrix))

# Create a data frame for ggplot
umap_df <- data.frame(
  X = umap_out$layout[, 1],
  Y = umap_out$layout[, 2],
  time_status = filtered_metadata$time_status
)

# Plot UMAP results with a solid white background
ggplot(umap_df, aes(x = X, y = Y, color = time_status)) +
  geom_point(size = 2) +
  theme_minimal() +
  ggtitle("UMAP Plot of Samples") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_line(color = "gray90")
  )

# Save the plot with a solid white background
ggsave(file.path(plots_dir, "umap_plot.png"))