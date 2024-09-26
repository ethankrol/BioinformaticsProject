#Load files if not already loaded
library(dplyr)
library(readr)
filtered_data_file <- "results/filtered_data.tsv"
filtered_metadata_file <- "results/filtered_metadata.tsv"

filtered_expression_df <- readr::read_tsv(filtered_data_file) %>%
  tibble::column_to_rownames("Symbol")

filtered_metadata <- readr::read_tsv(filtered_metadata_file)

# Filter out metadata where time_status is either one day or one month, and pull those sample IDs
# to use for indexing expression data

twentyfour_hour_samples <- filtered_metadata %>% dplyr::filter(time_status == "One_Day_After_Defeat") %>% dplyr::pull(refinebio_accession_code)

one_month_samples <- filtered_metadata %>% dplyr::filter(time_status == "One_Month_After_Defeat") %>% dplyr::pull(refinebio_accession_code)

# Use filtered metadata to match and set new dataframes with respective
# time_status values.

twentyfour_hour_expression <- filtered_expression_df[, twentyfour_hour_samples]

one_month_expression <- filtered_expression_df[, one_month_samples]

# Utilized this article for setting up wilcox.test and boxplot:
# https://library.virginia.edu/data/articles/the-wilcoxon-rank-sum-test

# Performing wilcox on each individual gene

wilcox_results_per_gene <- data.frame(gene = rownames(filtered_expression_df), p_value = NA)

# Loop through all the genes and perform wilcox rank-sum test on each gene

for (i in 1:nrow(filtered_expression_df)){
  gene_expression_twentyfour_hour <- as.numeric(twentyfour_hour_expression[i, ])
  gene_expression_one_month <- as.numeric(one_month_expression[i,])
  
  # Run wilcox rank-sum on the gene by getting a vector of all of the normalized
  # count values:
  
  wilcox_result <- wilcox.test(gene_expression_twentyfour_hour, gene_expression_one_month, exact = FALSE)
  
  wilcox_results_per_gene$p_value[i] <- wilcox_result$p.value
}

# Run clusterProfiler gene ontology
# Reference: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html

# BiocManager::install("clusterProfiler")

library(clusterProfiler)

significant_genes <- wilcox_results_per_gene$gene[wilcox_results_per_gene$p_value < 0.05]

gene_ontology_results <- enrichGO(gene = significant_genes,
                                  OrgDb = org.Mm.eg.db,
                                  keyType = "SYMBOL",
                                  ont = "ALL",
                                  pAdjustMethod = "BH",
                                  qvalueCutoff = .05)


# Take $result and convert to dataframe, and then take the top 20 results. Write
# both to files.

gene_ontology_results_df <- as.data.frame(gene_ontology_results)

readr::write_tsv(
  gene_ontology_results_df,
  file.path(results_dir, "wilcoxon_rank_sum_test_results.tsv")
)

top_enrichments <- gene_ontology_results_df %>%
  dplyr::arrange(pvalue) %>% 
  dplyr::slice(1:20)

readr::write_tsv(
  top_enrichments,
  file.path(results_dir, "wilcoxon_rank_sum_test_top_20_results.tsv")
)
