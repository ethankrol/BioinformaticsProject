chi_squared_test_cluster_sample <- c()

chi_10_vs_sample <- chisq.test(table(metadata_10_genes$cluster, filtered_metadata$time_status))
chi_squared_test_cluster_sample <- c(chi_squared_test_cluster_sample, chi_10_vs_sample$p.value)

# Chi-squared test for 100 genes
chi_100_vs_sample <- chisq.test(table(metadata_100_genes$cluster, filtered_metadata$time_status))
chi_squared_test_cluster_sample <- c(chi_squared_test_cluster_sample, chi_100_vs_sample$p.value)

# Chi-squared test for 1000 genes
chi_1000_vs_sample <- chisq.test(table(metadata_1000_genes$cluster, filtered_metadata$time_status))
chi_squared_test_cluster_sample <- c(chi_squared_test_cluster_sample, chi_1000_vs_sample$p.value)

# Chi-squared test for 5000 genes
chi_5000_vs_sample <- chisq.test(table(metadata_5000_genes$cluster, filtered_metadata$time_status))
chi_squared_test_cluster_sample <- c(chi_squared_test_cluster_sample, chi_5000_vs_sample$p.value)

# Chi-squared test for 10000 genes
chi_10000_vs_sample <- chisq.test(table(metadata_10000_genes$cluster, filtered_metadata$time_status))
chi_squared_test_cluster_sample <- c(chi_squared_test_cluster_sample, chi_10000_vs_sample$p.value)

adjusted_chi_squared_test_cluster_sample <- p.adjust(chi_squared_test_cluster_sample)

chi_squared_results_samples <- data.frame(
  Gene_Set = c("10 genes", "100 genes", "1000 genes", "5000 genes", "10000 genes"),
  Chi_Squared_Statistic = c(chi_10_vs_sample$statistic, chi_100_vs_sample$statistic, chi_1000_vs_sample$statistic, chi_5000_vs_sample$statistic, chi_10000_vs_sample$statistic),
  P_Value = chi_squared_test_cluster_sample,
  Adjusted_P_Value = adjusted_chi_squared_test_cluster_sample
)

write.csv(chi_squared_results_samples, "chi_squared_results_assignment1.csv", row.names = FALSE)
