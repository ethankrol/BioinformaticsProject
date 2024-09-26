source("config.r")

# Install TopGo Package
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("topGO")

# Install Libraries
library(topGO)
library(org.Mm.eg.db)

# Get all genes with name and padj
allVec <- as.vector(deseq_df$padj)
names(allVec) <- as.character(deseq_df$Gene)
# Filter out 5 NA, size from 22585 to 22580
length(allVec)
allVec <- na.omit(allVec)
length(allVec)

# Selection Function
selection <- function(allScore) { 
  return(allScore < 0.05) 
}

# Build topGO data
sampleGOdata <- new("topGOdata",
                    ontology = "BP", allGenes = allVec,
                    geneSel = selection, nodeSize = 10,
                    annot = annFUN.org, mapping = "org.Mm.eg.db",
                    ID = "symbol")

# Run Classic Enrichment Test
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

# List Significant GO Terms
# Use binary search to find the total GO Terms (6887)
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 6887)
print(allRes)

# write to tsv
write.table(allRes, 
            file = file.path(results_dir, "topGO_Result.tsv"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
