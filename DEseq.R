
# Exploratory data analysis


# Load libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
library(pheatmap)

# 1. Import count matrix

# Read counts (skip first line if it's a comment)
counts <- read.table("counts_clean.txt", header=TRUE, row.names=1, check.names=FALSE)

# Strip paths from column names
colnames(counts) <- sub(".*/(SRR[0-9]+)\\.sorted\\.bam", "\\1", colnames(counts))

# Now the column names are just sample IDs
colnames(counts)

dim(counts)
head(counts[,1:6])

# 2. Define sample metadata

sample_info <- data.frame(
  row.names = colnames(counts),
  genotype = c(
    "WT","WT","WT","WT","WT",
    "DKO","DKO","DKO","DKO",
    "WT","WT","WT",
    "DKO","DKO","DKO"
  ),
  infection = c(
    "Case","Case","Case","Case","Case",
    "Case","Case","Case","Case",
    "Control","Control","Control",
    "Control","Control","Control"
  )
)


# Verify order
all(rownames(sample_info) == colnames(counts))

# 3. Create DESeqDataSet

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ genotype + infection
)

# 4. Run DESeq normalization & modeling

dds <- DESeq(dds)

# 5. Variance-stabilizing transformation

vsd <- vst(dds, blind = TRUE)

# 6. PCA plot
# PCA plot
p <- plotPCA(vsd, intgroup = c("genotype", "infection"))  # ggplot object
print(p)

# Save PCA
ggsave("PCA_plot.png", p, width = 6, height = 5)

# PC1, PC2 values for each sample
pca_data <- plotPCA(vsd, intgroup = c("genotype", "infection"), returnData = TRUE)
pc1_table <- pca_data[, c("name", "PC1", "PC2")]
print(pc1_table)




# 7. Optional

# vsd is the variance-stabilized DESeqDataSet

# Perform PCA using prcomp directly on the assay
pca <- prcomp(t(assay(vsd)))  # transpose so samples are rows

# The loadings tells how much each gene contributes to each principal component
loadings <- pca$rotation  # genes × PCs

# Inspect top contributing genes for PC1
top_PC1_genes <- head(sort(abs(loadings[,1]), decreasing = TRUE), 20)
top_PC1_genes

#Differential expression analysis

# Example: DKO vs WT (ignoring infection)
res_genotype <- results(dds, contrast = c("genotype", "DKO", "WT"))
summary(res_genotype)

# Example: Case vs Control (within all samples)
res_infection <- results(dds, contrast = c("infection", "Case", "Control"))
summary(res_infection)

de_genes <- res_genotype[which(res_genotype$padj < 0.05), ]
nrow(de_genes)  # total DE genes
table(de_genes$log2FoldChange > 0)  # TRUE = up, FALSE = down

res_ordered <- res_genotype[order(res_genotype$padj), ]
res_ordered <- res_ordered[!is.na(res_ordered$padj), ]

top3_genes <- rownames(res_ordered)[1:3]
print(top3_genes)

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Extract counts for genes of interest
norm_counts[c("ENSMUSG00000032690", "ENSMUSG00000073409","ENSMUSG00000041827"), ]

