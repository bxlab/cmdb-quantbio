# Assignment 8: Bulk RNA-Seq
Assignment Date: Friday, Nov. 4, 2022 <br>
Due Date: Friday, Nov. 11, 2022 <br>

## Lecture

**Slides** are available here: Coming Soon


## Assignment Overview


## Data

The data we'll using today is from [this paper](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000590) that examined sex-specific gene expression in developing Drosophila embryos. The paper provides their own processed gene expression data, but we've re-processed the raw reads and generated our own transcript-level gene expression data, and this is what you'll be using for today's assignment.

The FPKM normalized gene expression data can be found here: `~/cmdb-quantbio/assignments/lab/bulk_RNA-seq/extra_data/dros_gene_expression.csv`. Copy this file to the `answers` directory you made for this assignment.

This file has FPKM data for over 34,000 transcripts of over 17,000 genes across 5 developmental stages (stages 10-14) for both male and female embryos. Take a look at the file with `less` to make sure you understand how the data is organized.


## Assignment

Create a python script for this assignment. Everything you'll need to do for this assignment will be done in this script.

### Step 0a: Read in the data

1. We recommend using `numpy` to work with this dataset (although you're welcome to use `pandas` if you feel comfortable doing so). First, you'll need read the `.csv` file into a `numpy` array. You can do so with the following code:

```
input_arr = np.genfromtxt("dros_gene_expression.csv", delimiter=',', names=True, dtype=None)
```

2. With this structured array, you have access to the transcript names (the rows) and the column names from the `.csv` file. Extract this info and store it in some variables. For the column names, you can use the following code. What can you do to only include the sample names in this?

```
col_names = [input_arr.dtype.names]
```

3. Then subset your input data to only include the FPKM values, excluding the transcript name info.

### Step 0b: Process input data

Before running any analyses, you'll need to process your data a little bit more:

1. Convert your structured 1D array into an unstructured 2D array that can be indexed with row and column indices. You should use the following code where `fpkm_values` is the array from Step0a:3, that only includes the FPKMS, not the transcript names:

```
import numpy.lib.recfunctions as rfn
fpkm_values_2d = rfn.structured_to_unstructured(fpkm_values, dtype=np.float)
```

You will want to use this `fpkm_values_2d` unstructured array to filter and transform your data in the next two substeps.

2. Subset your data to only the transcripts whose median expression is greater than 0. You can use `numpy.median()` to find the median expression of each transcript. This function has an `axis` argument that you need to set, which will determine whether you're correctly finding the median expression per transcript or instead finding the median expression per sample. [This page](https://stackoverflow.com/questions/22320534/how-does-the-axis-parameter-from-numpy-work) might help you determine what value you should assign to `axis`. After finding the median expression per transcript, you can use `numpy.where()` to subset your fpkm array appropriately. Make sure to also filter your transcript name variable.

3. Using your filtered array, apply a log2(FPKM + 0.1) transform to your data.


### Step 1: Clustering

We'd like to see if we can identify any broad patterns present in our gene expression data. To explore this, we're going to cluster the data, both by sample (i.e. identifying developmental timepoints whose gene expression patterns are similar across the genes we measured) as well as by transcript (i.e. identifying transcripts whose expression patterns are similar over the course of development).

To do this analysis, you'll be using the `dendrogram`, `linkage` and `leaves_list` functions from [`scipy`](http://docs.scipy.org/doc/scipy-0.14.0/reference/cluster.hierarchy.html). The documentation for SciPy isn't very helpful, but with some quick Googling you can find examples of how to use both of these tools.

1. Using `linkage` and `leaves_list`, cluster the filtered and log2 transformed gene expression data matrix (derived from `fpkm_values_2d` earlier) for both genes and samples based on their patterns of expression (so both the rows and columns of the matrix). You will find the [numpy transpose functionality](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.T.html) useful in order to cluster the columns.

2. Plot a heatmap of the clustered gene expression data.

3. Using the `dendrogram` function, create a dendrogram relating the samples to one another.


### Step 2: Differential Gene Expression

## Submission


## Additional Resources
