#set working directory. Make sure all the files are there
setwd("/Users/m255029/Desktop/9.07.23/miRNA_mRNA_Corelation_analysis_Aarti/DEseqanalysis/625_5p/test")
# Load necessary libraries
library(DESeq2)
library(ggplot2)

# Load count data 
countData <- read.csv("625_5p.csv", header = TRUE, row.names = 1)

# Load sample information 
sampleInfo <- read.csv("625_5pmetadata.csv", header = TRUE)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData, colData = sampleInfo, design = ~ Group)

#make sure that are same
ncol(countData)  # Number of columns in countData (samples)
nrow(sampleInfo) # Number of rows in sampleInfo (samples)

# Perform data normalization
dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized = TRUE)

# Log10 transformation. (Use anyone log10 or log2)
log10_normalized_counts <- log10(normalized_counts + 1)  # I am adding 1 to avoid log(0)
boxplot(log10_normalized_counts, main = "Log10 Transformed Counts Distribution hsa_miR_625_5p", xlab = "Samples", ylab = "Log10 Transformed Counts")

# Log2 transformation (Use anyone log10 or log2)
log2_normalized_counts <- log2(normalized_counts + 1)  # I am adding 1 to avoid log(0)
boxplot(log2_normalized_counts, main = "Log2 Transformed Counts Distribution hsa_miR_625_5p", xlab = "Samples", ylab = "Log2 Transformed Counts")


# Create a box plot of normalized counts for all samples
boxplot(normalized_counts, main = "Normalized Counts Distribution", xlab = "Samples", ylab = "Normalized Counts")


# Calculate the threshold for low expression (e.g., 5 counts)
threshold <- 5

# Calculate the percentage of samples with expression below the threshold for each gene
low_expr_percent <- rowMeans(counts(dds) < threshold)

# Filter out genes with low expression in more than 80% of samples
keep_genes <- low_expr_percent <= 0.80

# Subset the DESeqDataSet to keep only the genes meeting the criteria
dds_filtered <- dds[keep_genes,]

# Perform differential expression analysis (comparing GroupB to GroupA)
res <- results(dds_filtered, contrast = c("Group", "Up", "Down"))
result_file <- "deseq2_results.csv"
write.csv(res, result_file)
deseq2_results <- read.csv("deseq2_results.csv", header = TRUE, row.names = 1)

# Assuming you want to keep rows that exist in both countData and deseq2_results
common_rows <- intersect(rownames(countData), rownames(deseq2_results))
countData <- countData[common_rows, ]
deseq2_results <- deseq2_results[common_rows, ]
countDataWithResults <- cbind(countData, deseq2_results)

write.csv(countDataWithResults, "countData_with_deseq2_results.csv")

library(dplyr)

# Read the updated count data with DESeq2 results
countDataWithResults <- read.csv("countData_with_deseq2_results.csv", header = TRUE, row.names = 1)

# Set your significance threshold (e.g., adjusted p-value < 0.05)
alpha <- 0.05

# Filter the genes based on the adjusted p-value
significant_genes <- countDataWithResults %>%
  filter(padj < alpha)

# Save the significant genes to a new CSV file
write.csv(significant_genes, "significant_genes.csv")
# Filter for significant genes (optional)
alpha <- 0.05
sig_genes <- subset(res, padj < alpha)

# View the top differentially expressed genes
head(sig_genes)

write.csv(sig_genes, file.path(results_folder, "differential_expression_results.csv"))

# Save the updated count data with results to a new CSV file
write.csv(countDataWithResults, "countData_with_results.csv")

library(biomaRt)

# Specify the organism and dataset for Ensembl
organism <- "hsapiens"  # Replace with your organism of interest, e.g., "mmusculus" for mouse
ensembl_dataset <- "hsapiens_gene_ensembl"  # Replace with the appropriate dataset for your organism

# Connect to the Ensembl database
ensembl <- useEnsembl(biomart = "ensembl", dataset = ensembl_dataset)

# Read your DESeq analysis results CSV file (replace 'deseq_results.csv' with your file)
deseq_results <- read.csv("significant_genes.csv", header = TRUE)

# Extract the gene IDs from your DESeq results
gene_ids <- deseq_results$GeneID  # Replace 'GeneID' with the actual column name in your file

# Retrieve gene annotations using biomaRt
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "start_position", "end_position", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = ensembl
)

# Merge gene annotations with your DESeq results based on gene IDs
merged_results <- merge(deseq_results, gene_annotations, by.x = "GeneID", by.y = "ensembl_gene_id")

# Save the merged results to a new CSV file 
write.csv(merged_results, "merged_deseq_results.csv", row.names = FALSE)


#Script for Volcano plot
# Reset graphics parameters
par(mfrow=c(1,1))

# Make a basic volcano plot with a custom y-axis scale
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="hsa_miR_625_5p", xlim=c(-10,10), ylim=c(0, 6)))

# Add colored points: blue if log2FC < -1 or log2FC > 1, red if p-value < 0.05, and black if -1 <= log2FC <= 1
with(subset(res, log2FoldChange < -1 | log2FoldChange > 1), points(log2FoldChange, -log10(padj), pch=20, col="blue"))

# Filter and label significant genes with fold change > 1 or < -1
significant_genes <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))

# Highlight significant genes in red
points(significant_genes$log2FoldChange, -log10(significant_genes$padj), pch=20, col="red")
