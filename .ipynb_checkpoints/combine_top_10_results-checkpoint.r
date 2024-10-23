source("config.r")

library(dplyr)

gProfiler2_Gene_Ontology_top10_file <- file.path(results_dir, "gProfiler2_Gene_Ontology_top_10_gsea_go_enrichment_results.tsv")
clusterProfiler_top10_file <- file.path(results_dir, "SRP062829_clusterProfiler_HallmarkONTO_top_10_gsea_results.tsv")
topGO_MF_top10_file <- file.path(results_dir, "topGO_MF_top_10_results.tsv")
topGO_BP_top10_file <- file.path(results_dir, "topGO_BP_Top10_Result.tsv")
wilcoxon_top10_file <- file.path(results_dir, "wilcoxon_rank_sum_test_top_10_results.tsv")

gProfiler <- readr::read_tsv(gProfiler2_Gene_Ontology_top10_file) %>%
    dplyr::rename(id = term_id, p_value = p_value) %>%
    dplyr::mutate(method = "gProfiler2_Gene_Ontology") %>%
    dplyr::select(id, p_value, method) %>%
    dplyr::mutate(p_value = as.numeric(p_value))

clusterProfiler <- readr::read_tsv(clusterProfiler_top10_file) %>%
    dplyr::rename(id = ID, p_value = p.adjust) %>%
    dplyr::mutate(method = "clusterProfiler_Hallmark_Ontology") %>%
    dplyr::select(id, p_value, method) %>%
    dplyr::mutate(p_value = as.numeric(p_value))

topGO_MF <- readr::read_tsv(topGO_MF_top10_file) %>%
    dplyr::rename(id = GO.ID, p_value = classicFisher) %>%
    dplyr::mutate(method = "topGO_MF") %>%
    dplyr::select(id, p_value, method)%>%
    dplyr::mutate(p_value = gsub("[<>]", "", p_value)) %>%
    dplyr::mutate(p_value = as.numeric(p_value)) 


topGO_BP <- readr::read_tsv(topGO_BP_top10_file) %>%
    dplyr::rename(id = GO.ID, p_value = classicFisher) %>%
    dplyr::mutate(method = "topGO_BP") %>%
    dplyr::select(id, p_value, method)%>%
    dplyr::mutate(p_value = gsub("[<>]", "", p_value)) %>%
    dplyr::mutate(p_value = as.numeric(p_value))


wilcoxon <- readr::read_tsv(wilcoxon_top10_file) %>%
    dplyr::rename(id = ID, p_value = p.adjust) %>%
    dplyr::mutate(method = "wilcoxon_rank_sum_test_BP") %>%
    dplyr::select(id, p_value, method)%>%
    dplyr::mutate(p_value = as.numeric(p_value))


## combine results and sort by p value
all_results <- dplyr::bind_rows(gProfiler, clusterProfiler, topGO_MF, topGO_BP, wilcoxon) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(method = paste(method, collapse = ", "), min_p_value = min(p_value), resp_p_values = paste(p_value, collapse = ", ")) %>%
    dplyr::arrange(min_p_value)

## save the combined results
readr::write_tsv(all_results, file.path(results_dir, "combined_top_10_results.tsv"))