# Assignment 11: Single-Cell RNA-seq

Assignment Data: Friday, December 8, 2023<br />
Due Date: Friday, December 15, 2023

## Livecoding Script

[Livecoding Script](https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/scRNA-seq/slides_asynchronous_or_livecoding_resources/livecoding.py)

## Mamba Environment

We will use `scanpy` for this lab, a fairly comprehensive Python package for scRNAseq analysis. We’ll create a mamba environment that we can use for this lab called `scanpy`. Use the command:

```bash
mamba create -n scanpy scanpy openblas "matplotlib<3.7" leidenalg "pandas<2.0.1" -y
```

Then you can activate the mamba environment with the command `mamba activate scanpy`.

## Getting the Data

We will be looking at a dataset containing ~3,000 peripheral blood mononuclear cells (PBMCs) from a healthy human donor. This was produced using the 10x technology and provided by 10X Genomics.

**What’s already been done:** the CellRanger package from 10x was used to align and count the reads. Reads were de-duplicated using UMIs and separated into “cells” based on barcodes, and then aligned to a transcriptome (using the STAR aligner). The result is a cell x gene matrix, which has been stored in a sparse matrix format.

Once you’re in the directory where you want the data, you can use the following command to download it:

```bash
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

This will create a folder, `filtered_gene_bc_matrices/hg19/` with 3 files:

- barcodes.tsv
- genes.tsv
- matrix.mtx

You will also need the livecoding excercise code. You can download the code [here](https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/scRNA-seq/slides_asynchronous_or_livecoding_resources/livecoding.py). Run it to produce the 2 input files you will need for the homework assignment, `filtered_data.h5` and `variable_data.h5`.

## The Assignment

### Getting the data into Scanpy

To get you started, all access to Scanpy is typically through a module called `scanpy` which is imported under the name `sc` for convenience. We then load the count matrix into a table, which is [an instance of `AnnData`](https://scanpy.readthedocs.io/en/latest/usage-principles.html#anndata).

```python
import scanpy as sc
# Read the 10x dataset filtered down to just the highly-variable genes
adata = sc.read_h5ad("variable_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 
```

For the rest of this lab, refer to the [Scanpy documentation](https://scanpy.readthedocs.io/en/stable/), and specifically the [API documentation](https://scanpy.readthedocs.io/en/stable/api.html).

### Step 1: Clustering

You will need to calculate a neighborhood graph to determine the neighborhood space for each cell. This is a pre-processing function and you should use the parameters `n_neighbors=10` and `n_pcs=40`.

Next, use `leiden` clustering to identify clusters in the data. Produce t-SNE and UMAP plots showing the clusters. Note that to create the UMAP transformation, you need to specify `maxiter`. I suggest 900.

### Step 2: Distinguishing Genes

Identify and plot genes that distinguish each cluster. Use both the Wilcoxon and logistic regression approaches, implemented through the `rank_genes_groups` function. (Again, see the `sc.tl` module to actually perform the ranking and `sc.pl` for plotting). You will need to use the arguments `groupby='leiden'` and `use_raw=True` for both rankings.

### Step 3: Reload Missing Genes

You may receive errors trying to plot the data for certain genes that you pull from the highly-ranked genes. That's because the ranking is done across all genes, Because some genes are filtered out at the "highly-variable gene" filtering step, despite showing cell-type specificity, you should load the dataset `filtered_data.h5` which has the genes that were removed by the highly-variable gene filtering. In order to keep the leiden clusters and UMAP or t-SNE mappings, you will need to copy them from your clustered and transformed data. To do this, you can copy the values to a variable and then copy that variable to the `filtered_data.h5` dataset.

```python
leiden = adata.obs['leiden']
umap = adata.obsm['X_umap']
tsne = adata.obsm['X_tsne']
adata = sc.read_h5ad('filtered_data.h5')
adata.obs['leiden'] = leiden
adata.obsm['X_umap'] = umap
adata.obsm['X_tsne'] = tsne
```

You can do this as a single script or by creating 2 scripts and saving the data at the end of the first script and loading it in the second, to avoid having to rerun the clustering steps.

### Step 4: Cell Types?

Now the fun part.

Using your knowledge (or the knowledge of hematopoeisis aficionados in your cohort, or Google), identify some marker genes that should distinguish different blood cell types. You must identify at least 4 cell types. There are many resources online for identifying relationships between marker genes and cell types. Take a look at [this database](http://betsholtzlab.org/VascularSingleCells/database.html) for an example of ways to identify cell types. You may need multiple markers to distinguish between clusters.

1. Support plots that provide evidence for your cell type assignments/what you used to diagnose/decide on cell types.
  * You can color UMAP and t-SNE plots by any gene of your choice, which is helpful for visualizing which clusters are enriched for which genes, and which clusters might correspond to a specific blood cell type.
  * Alternatively, you can also produce `dotplots` and `clustermaps` that allow you to see how a specific set of genes are associated with your clusters. Also, stacked violin plots, etc…
2. Besides these support plots, make an overall t-SNE or UMAP plot that labels your clusters with the cell types you think they mostly represent, either on the plot or in a legend. Make sure to provide the support plots you made in order to establish your labeling. See [this tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/plotting/core.html) for an example of how to apply labels.

## Submit

- Python code
- t-SNE and UMAP plots of clusters from step 1
- Plots for genes that distinguish clusters (t-test and logistic regression) from step 2
- Support plots that provide evidence for your cell type assignments (Step 4.1)
- The overall t-SNE or UMAP plot with at least 4 cell types labeled (Step 4.2)

## Advanced Excercise

Identify all clusters (there should be 8). Create a dotplot with a gene as specific to each cluster as you can find (8 clusters, 8 genes). Create a multi-panel UMAP or t-SNE plot colored by the expression of each gene as well as the leiden clusters.