# Set working directory
setwd("C:/Users/m215200/Desktop/")

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("DESeq2", "pheatmap", "EnhancedVolcano"), ask = FALSE)
install.packages("ggplot2")

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)

# Load count data
countData <- read.csv("625_5p.csv", header = TRUE, row.names = 1)

# Load sample metadata
sampleInfo <- read.csv("625_5pmetadata.csv", header = TRUE, row.names = 1)

# Ensure column names match between countData and sampleInfo
countData <- countData[, rownames(sampleInfo)]

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleInfo,
                              design = ~ condition)  # Replace 'condition' with your actual group column

# Filter out low-count genes (optional but recommended)
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Get results (change contrast if needed)
res <- results(dds)

# Order by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Save results to CSV
write.csv(as.data.frame(resOrdered), file = "deseq2_results.csv")

# Log-normalize counts for visualization
rld <- rlog(dds, blind = FALSE)

# PCA plot
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: RNA-Seq Samples (DESeq2 rlog)")

ggsave("pca_plot.png", plot = p, width = 6, height = 4)

# MA plot
png("ma_plot.png", width = 700, height = 500)
plotMA(res, ylim = c(-5, 5), main = "MA Plot")
dev.off()

# Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                caption = "Log2FC cutoff = 1.5, p < 0.05")

ggsave("volcano_plot.png", width = 7, height = 6)
