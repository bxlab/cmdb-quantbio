# install new packages
BiocManager::install("DESeq2")
BiocManager::install("vsn")

# Load libraries we'll need
library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(ggfortify)

# Load tab-separated data file
data = readr::read_tsv('salmon.merged.gene_counts.tsv')

# Change gene names into row names
data = column_to_rownames(data, var="gene_name")

# Get rid of gene id column
data = data %>% select(-gene_id)

# Change data to integers
data = data %>% mutate_if(is.numeric, as.integer)

# Remove low coverage samples
data = data[rowSums(data) > 100,]

# Pull out broad region samples
broad = data %>% select("A1-3_Rep1":"P1-4_Rep3")

# Create metadata tibble with tissues and replicate numbers based on sample names
broad_metadata = tibble(tissue=as.factor(c("A1_3", "A1_3", "A1_3", "Cu_LFC_Fe", "Cu_LFC_Fe", "Cu_LFC_Fe", "P1_4", "P1_4", "P1_4")),
                        rep=as.factor(c(1, 2, 3, 1, 2, 3, 1, 2, 3)))

# Create a DESeq data object
broaddata = DESeqDataSetFromMatrix(countData=as.matrix(broad), colData=broad_metadata, design=~tissue)

# Plot variance by average
meanSdPlot(assay(broadLogdata))

# Log transform data
broadLogdata = normTransform(broaddata) # log(1 + data)

# Plot log-transformed data variance by average
meanSdPlot(assay(broadLogdata))

# Create PCA data
broadPcaData = plotPCA(broadLogdata,intgroup=c("rep","tissue"), returnData=TRUE)

# Plot PCA data
ggplot(broadPcaData, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5)

# Batch-correct data to remove excess variance with variance stabilizing transformation
broadVstdata = vst(broaddata)

# Plot variance by average to verify the removal of batch-effects
meanSdPlot(assay(broadVstdata))

# Perform PCA and plot to check batch-correction
broadPcaData = plotPCA(broadVstdata,intgroup=c("rep","tissue"), returnData=TRUE)
ggplot(broadPcaData, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5)

# Convert into a matrix
broadVstdata = as.matrix(assay(broadVstdata[sds>1,]))

# Find replicate means
combined = broadVstdata[,seq(1, 21, 3)]
combined = combined + broadVstdata[,seq(2, 21, 3)]
combined = combined + broadVstdata[,seq(3, 21, 3)]
combined = combined / 3

# Use replicate means to filter low variance genes out
filt = rowSds(combined) > 1
broadVstdata = broadVstdata[filt,]

# Plot expression values with hierarchical clustering
heatmap(broadVstdata, Colv=NA)

# Perform new hierarchical clustering with different clustering method
distance = dist(broadVstdata)
Z = hclust(distance, method='ave')

# Plot expression values with new hierarchical clustering
heatmap(broadVstdata, Colv=NA, Rowv=as.dendrogram(Z))

# Set seed so this clustering is reproducible
set.seed(42)

# Cluster genes using k-means clustering
k=kmeans(broadVstdata, centers=12)$cluster

# Find ordering of samples to put them into their clusters
ordering = order(k)

# Reorder genes
k = k[ordering]

# Plot heatmap of expressions and clusters
heatmap(broadVstdata[ordering,],Rowv=NA,Colv=NA,RowSideColors = RColorBrewer::brewer.pal(12,"Paired")[k])

# Save heatmap
png("heatmap.jpg")
heatmap(broadVstdata[ordering,],Rowv=NA,Colv=NA,RowSideColors = RColorBrewer::brewer.pal(12,"Paired")[k])
dev.off()

# Pull out gene names from a specific cluster
genes = rownames(broadVstdata[k == 9,])

# Same gene names to a text file
write.table(genes, "cluster_genes.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

