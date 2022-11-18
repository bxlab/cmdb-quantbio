# Assignment 10: Single-Cell RNA-seq

Assignment Data: Friday, November 18, 2022<br />
Due Date: Friday, December 2, 2022

## Lecture

[Lecture slides](https://docs.google.com/presentation/d/1n6T5zQXMepT8sP4ua8CS_hiGflxIm68DTZ0SczzJuHI/edit?usp=sharing)

## Conda Environment

We will use `scanpy` for this lab, a fairly comprehensive Python package for scRNAseq analysis. We’ll create a conda environment that we can use for this lab called `scanpy`. Use the command `conda create -n scanpy scanpy openblas matplotlib leidenalg multicore-tsne -y`. Then you can activate the conda environment with the command `conda activate scanpy`.

## Getting the Data

We will be looking at a dataset containing ~10,000 brain cells from an E18 mouse. This was produced using the 10x technology, using their most recent (v3) chemistry.

**What’s already been done:** the CellRanger package from 10x was used to align and count the reads. Reads were de-duplicated using UMIs and separated into “cells” based on barcodes, and then aligned to a transcriptome (using the STAR aligner). The result is a cell x gene matrix, which has been stored in a binary format (an hdf5 file).

Once you’re in the directory where you want the data, you can use the following command to download it:

```bash
curl https://bx.bio.jhu.edu/data/msauria/cmdb-lab/neuron_10k_v3_filtered_feature_bc_matrix.h5 \
--output neuron_10k_v3_filtered_feature_bc_matrix.h5
```

## The Assignment

### Getting the data into Scanpy

To get you started, all access to Scanpy is typically through a module call `scanpy` which is imported under the name `sc` for convenience. We then load the count matrix into a table, which is an instance of `AnnData`.

```python
import scanpy as sc
# Read 10x dataset
adata = sc.read_10x_h5("neuron_10k_v3_filtered_feature_bc_matrix.h5")
# Make variable names (in this case the genes) unique
adata.var_names_make_unique()
```

For the rest of this lab, refer to the [Scanpy documentation](https://scanpy.readthedocs.io/en/stable/), and specifically the [API documentation](https://scanpy.readthedocs.io/en/stable/api.html).

### Step 1: Filtering

Filtering tools are largely under the `sc.pp` module. We suggest using the *Zheng et al. 2017* filtering approach. Produce a PCA plot before and after filtering (see the `sc.tl` module to actually perform the PCA and `sc.pl` for plotting).

### Step 2: Clustering

Use `leiden` clustering to identify clusters in the data. Produce t-SNE and UMAP plots showing the clusters. Note that to create the UMAP transformation, you need to specify `maxiter`. I suggest 1000.

### Step 3: Distinguishing Genes

Identify and plot genes that distinguish each cluster. Use both the t-test and logistic regression approaches, implemented through the `rank_genes_groups` function.

### Step 4: Cell Types?

Now the fun part.

Using your knowledge, identify some marker genes that should distinguish different brain cell types. You must identify at least 6 cell types. There are many resources online for identifying relationships between marker genes and cell types. Take a look at [this database](http://betsholtzlab.org/VascularSingleCells/database.html) for an example of ways to identify cell types.

1. You can color UMAP and t-SNE plots by any gene of your choice, which is helpful for visualizing which clusters are enriched for which genes, and which clusters might correspond to a specific brain cell type.
2. Alternatively, you can also produce `dotplots` and `clustermaps` that allow you to see how a specific set of genes are associated with your clusters. Also, stacked violin plots, etc…

Besides these support plots, make an overall t-SNE or UMAP plot that labels your clusters with the cell types you think they mostly represent. Make sure to provide the support plots you made in order to establish your labeling. See [this tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/plotting/core.html) for an example of how to apply labels.

## Submit

- Python code
- PCA plots before and after filtering from step 1
- t-SNE and UMAP plots of clusters from step 2
- Plots for genes that distinguish clusters (t-test and logistic regression) from step 3
- The overall t-SNE or UMAP plot with at least 6 cell types labeled
- Support plots that provide evidence for your cell type assignments
