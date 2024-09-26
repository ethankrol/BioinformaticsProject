source("config.r")

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("topGO")

library(topGO)
library(org.Mm.eg.db)

allVec <- as.vector(deseq_df$padj)
names(allVec) <- as.character(deseq_df$Gene)
length(allVec)
allVec <- na.omit(allVec)
length(allVec)

selection <- function(allScore) { 
  return(allScore < 0.05) 
}

sampleGOdata <- new("topGOdata",
                    ontology = "MF", allGenes = allVec,
                    geneSel = selection, nodeSize = 10,
                    annot = annFUN.org, mapping = "org.Mm.eg.db",
                    ID = "symbol")

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 1346)
print(allRes)

write.table(allRes, 
            file = file.path(results_dir, "topGO_Result_MF.tsv"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
