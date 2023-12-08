# Single-Cell RNA-seq

## Livecoding Script

[Livecoding Script](https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/scRNA-seq/slides_asynchronous_or_livecoding_resources/livecoding.py)

## Mamba Environment

We will use `scanpy` for this lab, a fairly comprehensive Python package for scRNA-seq analysis. Create a mamba environment that you will use for this lab called `scanpy` using the following command:

```bash
mamba create -n scanpy scanpy openblas "matplotlib<3.7" leidenalg "pandas<2.0.1" -y
```

Then, activate the mamba environment with the command `mamba activate scanpy`.<br><br>

## Scanpy documentation

For the rest of this lab, refer to the [Scanpy documentation](https://scanpy.readthedocs.io/en/stable/), and specifically the [API documentation](https://scanpy.readthedocs.io/en/stable/api.html).

This documentation is broken into several sections, grouping functions by stage of analysis. This grouping is reflected in the prefix of each function:

- Preprocessing functions (filtering, qc, normalizing, etc.) begin with `pp.`
- Tool functions (analysis, embedding, clustering, etc.) begin with `tl.`
- Plotting functions (visualizing t-SNE, UMAP, PCA, expression, ranking, etc.) begin with `pl.`
- Reading and writing functions have no prefix but are called directly from the `scanpy` package

## Getting the Data

You will be looking at a dataset comprising scRNA-seq data of ~3,000 peripheral blood mononuclear cells (PBMCs) from a healthy human donor. This data was produced by 10X Genomics using their dedicated tools.

**What’s already been done:** the CellRanger package from 10x was used to align and count the reads. Reads were de-duplicated using UMIs (unique molecular identifiers) and separated into “cells” based on barcodes, and then aligned to a transcriptome (using the STAR aligner). The result is a cell x gene matrix of expression, which has been stored in a sparse matrix format.

Once you’re in the directory where you want the data, you can use the following command to download it:

```bash
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

This will create a folder, `filtered_gene_bc_matrices/hg19/` with 3 files:

- `barcodes.tsv`
- `genes.tsv`
- `matrix.mtx`

### Live-coding preprocessing

You will also need the live-coding excercise code. You can download the complete code [here](https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/scRNA-seq/slides_asynchronous_or_livecoding_resources/livecoding.py).

The live-coding script loads the pre-aligned and deduplicated reads into scanpy. It then preforms the following steps:

1. Filters the reads for under-represented genes and cells
2. Removes mitochondrial genes
3. Normalizes and log transforms expression levels using a pseudo-count of one
4. Saves this version of the data under the name `filtered_data.h5`
5. Remove genes that are not considered 'highly variable'
6. Regresses out the influence of the total transcript count per gene and the percent of mitochondrial gene expression
7. Renormalized expression levels with a maximum expression value of 10
8. Calculates principal components of the data
9. Saves the fully-filtered version of the data under the name `variable_data.h5`

Using the three files above, run the live-coding script to produce the 2 input files you will need for the homework assignment: `filtered_data.h5` and `variable_data.h5`.<br><br>

## Exercises

### Exercise 0: Getting the data into Scanpy

You will first need to load the counts matrix from the live-coding exercise into a special kind of table, which is [an instance of Scanpy's `AnnData` class](https://scanpy.readthedocs.io/en/latest/usage-principles.html#anndata).

You will be modifying this table throughout the homework as you run different analyses, and using it for plotting, all using scanpy functions.

Load the counts matrix into a scanpy `AnnData` object using the code below:

```python
import scanpy as sc
# Read the 10x dataset filtered down to just the highly-variable genes
adata = sc.read_h5ad("variable_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 
```

<span style="color: red">**NOTE:**</span> Most `scanpy` functions will modify your `adata` object *in place*. This means you generally don't need to create new variables when you run `scanpy` functions. You can simply run the function on the `adata` object, and it will modify it/add data to it as necessary.<br><br>

### Exercise 1: Clustering

A good place to start with scRNA-seq data is to try to cluster cells based on their expression profiles. In this exercise you will be doing just that: clustering cells based on their expression profiles, and then visualizing those clusters, using tSNE and UMAP.<br><br>

#### **Step 1.1**: Computing a neighborhood graph

Before running clustering, scanpy needs to construct what is called a "neighborhood" graph, which essentially records each cell's "nearest neighbors": the cells whose expression looks the most similar to the focal cell.

Using the [scanpy documentation](https://scanpy.readthedocs.io/en/stable/api.html), find a **pre-processing** function that computes a neighborhood graph from `adata`. Run the function you found using the parameters `n_neighbors=10` and `n_pcs=40`.<br><br>

#### **Step 1.2**: Leiden clustering

Next, you will use leiden clustering (this is just a specific clustering algorithm, there are many others) to identify clusters in the data. This clustering algorithm will use the neighborhood graph you computed in Step 1.1 to actually produce the clusters.

As before, use the [scanpy documentation](https://scanpy.readthedocs.io/en/stable/api.html) to find a **tool** function that will run leiden clustering on your data. Run the function you found; you shouldn't need any additional parameters.<br><br>

#### **Step 1.3**: Visualizing clusters

In this step, you will be visualizing the clusters you identified in the previous steps using two separate approaches: UMAP and t-SNE. UMAP and t-SNE are both dimensionality reduction tools for visualizing structure in your data, just like PCA, but they make fewer assumptions about how the data actually is structured.

Before you can produce UMAP and tSNE plots, you'll actually need to run the UMAP and tSNE algorithms, which essentially "embed" your high-dimensional data into two dimensions.

Use the [scanpy documentation](https://scanpy.readthedocs.io/en/stable/api.html) to find a **tool** function to run the UMAP algorithm on your data. Note that to create the UMAP transformation, you need to specify `maxiter`. We suggest 900.

Now, find a **tool** function to run the tSNE algorithm on your data. You should not need any additional parameters.

Now that you've run the UMAP and t-SNE algorithms, you'll want to plot the output. Helpfully, scanpy's plotting functions actually work with Matplotlib's `fig, ax = plt.subplots()` functionality: when you run a scanpy plotting function, you can use the `ax` argument to plot to a specific matplotlib `axes` object. Create a figure showing the UMAP and t-SNE visualizations of your data as follows:
1. Create a two panel figure using `fig, axes = plt.subplots(ncols=2)`
2. Find a scanpy **plotting** function to plot your UMAP results. You'll want to specify the `color` argument to color the cells by their leiden cluster (`color='leiden'`), the `ax` argument to choose where to actually plot the data (`ax = axes[0]`), and the `title` argument to set the panel title. You'll also want to set `show=False`.
3. Find a scanpy **plotting** function to plot your t-SNE results. You'll want to specify the `color` argument to color the cells by their leiden cluster (`color='leiden'`), the `ax` argument to choose where to actually plot the data (`ax = axes[1]`), and the `title` argument to set the panel title. You'll also want to set `show=False`.

Label your figure appropriately and save it to a file. You will be uploading it with the assignment.<br><br>

### Exercise 2: Identifying cluster marker genes

Now that you've identified clusters in your data, ideally you'd like to be able to determine what cell types those clusters correspond to. The first step to doing that is to identify the "marker genes" that seem to define each cluster: which genes are specifically upregulated in each cluster, relative to the other clusters. In this exercise, you'll be identifying each cluster's marker genes.<br><br>

#### **Step 2.1**: Ranking genes in each cluster

Scanpy has a great **tool** function for identifying and ranking potential marker genes in each cluster: `rank_genes_groups()`. This function has multiple different methods for ranking genes; you'll be using and comparing two: 1) Wilcoxon rank-sum and 2) logistic regression.

First: using the `sc.tl.rank_genes_groups()` function, rank marker genes in your data using the **Wilcoxon rank-sum** method. You'll want to specify the `method` argument to use the Wilcoxon rank-sum method, as well as the `groupby` argument to use the leiden clusters you identified earlier (`groupby='leiden'`). You'll also want to set `use_raw=True` and `copy=True`. **Store the results of running this function in a new `wilcoxon_adata` variable.** This will essentially be a copy of the `adata` object, but with the Wilcoxon `rank_genes_groups` data added.

Second: using the `sc.tl.rank_genes_groups()` function, rank marker genes in your data using the **logistic regression** method. You'll want to specify the `method` argument to use the logistic regression method, as well as the `groupby` argument to use the leiden clusters you identified earlier (`groupby='leiden'`). You'll also want to set `use_raw=True` and `copy=True`. **Store the results of running this function in a new `logreg_adata` variable.** This will essentially be a copy of the `adata` object, but with the Logistic regression `rank_genes_groups` data added.<br><br>

#### **Step 2.2**: Visualizing marker genes

Helpfully, scanpy has another **plotting** function for visualizing the gene rankings from `sc.tl.rank_genes_groups()`.

First, using the `sc.pl.rank_genes_groups()` function with the `wilcoxon_adata` object you created in the previous step, plot the top 25 genes for each cluster when using the **Wilcoxon rank-sum method**. As before, you can use Matplotlib's `fig, ax = plt.subplots()` functionality to make your life easier. You'll want to specify the `n_genes`, `ax` and `title` arguments. You'll also want to set `sharey=False`, `show=False`, and `use_raw=True`. Label your figure appropriately and save it to a file. You will be uploading it with the assignment.

Second, using the `sc.pl.rank_genes_groups()` function with the `logreg` object you created in the previous step, plot the top 25 genes for each cluster when using the **logistic regression method**. As before, you can use Matplotlib's `fig, ax = plt.subplots()` functionality to make your life easier. You'll want to specify the `n_genes`, `ax` and `title` arguments. You'll also want to set `sharey=False`, `show=False`, and `use_raw=True`. Label your figure appropriately and save it to a file. You will be uploading it with the assignment.<br><br>

### Exercise 3: Identifying cluster cell types

Now that you've identified marker genes for each cluster in your data set, you want to assign cell-types to each cluster, based on those marker genes.<br><br>

#### **Step 3.1**: Reload Missing Genes

In the next couple of steps, you may receive errors trying to plot the data for certain genes that you pull from the highly-ranked genes. That's because the ranking is done across all genes. Because some genes are filtered out at the "highly-variable gene" filtering step, despite showing cell-type specificity, you should load the dataset `filtered_data.h5` which has the genes that were removed by the highly-variable gene filtering. In order to keep the leiden clusters and UMAP or t-SNE mappings, you will need to copy them from your clustered and transformed data. To do this, you can copy the values to a variable and then copy that variable to the `filtered_data.h5` dataset.

```python
leiden = adata.obs['leiden']
umap = adata.obsm['X_umap']
tsne = adata.obsm['X_tsne']
adata = sc.read_h5ad('filtered_data.h5')
adata.obs['leiden'] = leiden
adata.obsm['X_umap'] = umap
adata.obsm['X_tsne'] = tsne
```

You can do this as a single script or by creating 2 scripts and saving the data at the end of the first script and loading it in the second, to avoid having to rerun the clustering steps, like this:

```python
adata.write('filtered_clustered_data.h5')
```

and load in a new script:

```python
adata = sc.read_h5ad("filtered_clustered_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 
```
<br>

#### **Step 3.2**: Matching genes to cell types

Now the fun part.

Identify some marker genes in your data set that should distinguish different blood cell types. Your *goal* is to try to identify at least 3 distinct cell types in your data based on the marker genes from exercise 2. You can do this using your knowledge, the knowledge of hematopoeisis aficionados in your cohort, or the power of the internet. There are many resources online for identifying relationships between marker genes and cell types. Take a look at [this database](http://betsholtzlab.org/VascularSingleCells/database.html) for an example of ways to identify cell types. You may need multiple marker genes to distinguish between clusters.

Once you have a short list of genes that should distinguish cell-types in your data, produce ONE figure that supports your hypothesis that these genes should be marking specific clusters in your dataset. You have several options (and [this link](https://scanpy-tutorials.readthedocs.io/en/latest/plotting/core.html) may help you):
1. A set of UMAP or t-SNE plots, colored by expression of a specific gene. You can color the plots by a gene by using the `name` argument. Create a multi-panel figure (one panel per marker gene) showing the expression of these genes in either a UMAP or t-SNE plot. Compare your results to your plot from step 1.3 to validate your hypothesis (no need to write anything down).
2. A `dotplot` showing expression of your chosen marker genes in each cluster
3. A `clustermap` showing expression of your chosen marker genes in each cluster
4. Other (e.g. stacked violin plots)

Choose one approach and produce the corresponding figure. Make sure that the marker genes you chose do indeed match your expectations. Label your figure appropriately and save it to a file. You will be uploading it with the assignment.<br><br>

#### **Step 3.3**: Putting it all together

Now that you've confirmed that the marker genes you've chosen do indeed "mark" the three clusters you expected them to, you want to update your cluster labels with the cell-types corresponding to those clusters.

Using `adata.rename_categories()` rename the three clusters you identified in your `adata` object to the cell types you think they are (based on the marker genes).

Finally, make an overall t-SNE or UMAP plot that labels your clusters with the cell types you think they mostly represent, either on the plot or in a legend. Label your figure appropriately and save it to a file. You will be uploading it with the assignment.<br><br>

## Submission

1. Python script with all code for the assignment (**2 points total**)
  * Compute neighborhood graph (**0.5 point**)
  * Cluster data (**0.5 point**)
  * Rank marker genes using Wilcoxon rank-sum and logistic regression methods (**1 point**)
2. Pretty plots (**8 points total**)
  * UMAP/t-SNE multi-panel plot from step 1.3 (**2 points**)
  * Wilcoxon-based rank_genes_groups plot from step 2.2 (**1 point**)
  * Logreg-based rank_genes_groups plot from step 2.2 (**1 point**)
  * Chosen support plot from step 3.2 (**2 point**)
  * Final labeled UMAP or t-SNE plot from step 3.3 (**2 point**)

**Total Points: 10**<br><br>

## Advanced Excercise

Identify *all* clusters (there should be 8), in addition to the three you chose in step 3.2. Create a dotplot with a gene as specific to each cluster as you can find (8 clusters, 8 genes). Create a multi-panel UMAP or t-SNE plot colored by the expression of each gene as well as the leiden clusters.<br><br>
