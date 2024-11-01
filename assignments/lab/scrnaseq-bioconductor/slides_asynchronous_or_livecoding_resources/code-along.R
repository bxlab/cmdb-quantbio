#
# Basic single cell analysis inspired by Bioconductor OSCA Quick start 
# - https://bioconductor.org/books/3.19/OSCA.intro/analysis-overview.html#quick-start-simple
#



### Load packages

library( "DropletUtils" )  # Utilities for handling 10x Genomics data
library( "scater" )        # Tools for quality control and visualization
library( "scran" )         # Methods for single cell data analysis



### Load data

# wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
sce <- read10xCounts( "filtered_gene_bc_matrices/hg19/" )
sce



### Explore SingleCellExperiment

# cell metadata
colData(sce)
# convert from DataFrame to data.frame
as.data.frame(colData(sce))

# gene metadata
rowData(sce)
# change from ID to Symbol
rownames(sce) <- rowData(sce)$Symbol
sce

# subset matrix and SingleCellExperiment with [row,col] 
assay(sce)[40:50,1:10]
assay(sce[40:50,1:10])
assay(sce[40:50,])[,1:10]

# sum up rows or columns for insight
genecounts <- rowSums(assay(sce))
summary( genecounts )
head( sort(genecounts, decreasing=TRUE ) )



### Basic processing

# Library size normalization and log-transformation
sce <- logNormCounts( sce )

# Quantify per-gene variation
dec <- modelGeneVar(sce)
# Select highly variable genes
hvg <- getTopHVGs(dec, prop=0.1)

# Principal component analysis
set.seed(1234) # Set seed for reproducibility
sce <- runPCA( sce, subset_row=hvg )
plotReducedDim( sce, "PCA" )



### Clustering

# Principal components are used for clustering
colLabels(sce) <- clusterCells(sce, use.dimred='PCA')
as.data.frame(colData(sce))
plotReducedDim( sce, "PCA", colour_by="label" )



### Additional visualizations

# TSNE and UMAP primarily for visualization
set.seed(1234) # Set seed for reproducibility
sce <- runTSNE(sce, dimred = 'PCA')
plotReducedDim( sce, "TSNE", colour_by="label" )

# Plotting distributions is a good complement
rownames(sce) <- rowData(sce)$Symbol
goi <- c( "MS4A1", "CD14", "CD8A", "IL7R" )
plotExpression( sce, features=goi, x="label" )



### Finding marker genes

# All-vs-all comparison of clusters
marker.info <- scoreMarkers(sce, colData(sce)$label)
marker.info

# Sort based on AUC
chosen <- marker.info[["11"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])

plotReducedDim( sce, "TSNE", colour_by="CD79A" )
plotExpression( sce, "CD79A", x="label" )
