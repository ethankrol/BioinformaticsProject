library(pheatmap)

scaled_top_5000_genes <- t(scale(t(top_5000_gene_data)))

filtered_metadata <- readr::read_tsv("/results/filtered_metadata.tsv")

predictive_data <- readr::read_tsv("results/predictive_model_labels_5000_genes.tsv")
predictive_data$time_status <- predictive_data$RFC_label

heatmap_annotations <- data.frame(
  SVM = as.factor(predictive_data$SVM_label),
  RFC = as.factor(predictive_data$RFC_label),
  LogisticRegression = as.factor(predictive_data$LogisticRegression_label),
  KNN = as.factor(predictive_data$KNN_label),
  NaiveBayes = as.factor(predictive_data$NaiveBayes_label),
  SampleGroups= as.factor(predictive_data$time_status)
)

rownames(heatmap_annotations) <- colnames(scaled_top_5000_genes)

heatmap_plot <- pheatmap(
  scaled_top_5000_genes,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_col = heatmap_annotations,
  main = "Heatmap of Predictive Modeling Genes and their Correlated Models"
)
