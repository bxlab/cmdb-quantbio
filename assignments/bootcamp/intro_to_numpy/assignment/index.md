## Using Numpy matrices and indexing to find marker genes in GTEx data

For a description of GTEx, see the [GTEx Portal](https://gtexportal.org/home/aboutAdultGtex).

You will be using a summary dataset that has the median expression levels across samples for each gene and tissue, named `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct`.

Please submit your answers as a `.py` script with comments explaining the logic of each code block, a `.R` script for plotting the results, a `.pdf` plot (step #6), and a text file with the results from step #10. (Remember that comment lines start with `#` and are ignored by Python and R). For questions regarding interpretation or discussion of your results, please include your answers as comments interspersed with your code.

1. Load the GTEx file using the standard file `open` function discussed in previous sections.
- In order to load the data, first look at the file formatting using `less GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct` on the command line. Note that there are lines at the beginning that you want to skip. You can use the filestream method `readline()` to read in one line at a time. This will allow you to read in header lines without saving the results.
- There are four types of information you will want to capture: tissue names, gene IDs, gene names, and expression levels. Each field is separated by a `\t`, so make sure to specify that when using the `split()` method.
- You should end up with a list of tissue names, a list of gene IDs, a list of gene names, and a nested list of expression levels

2.  Convert each of your data lists into a numpy array. In the case of text names, numpy will determine what datatype is appropriate. For your expression data, you will need to specify the data type as `float` using the argument `dtype=float`.

Why do you need to tell numpy what type of data the expression data are but not for the other data lists?

3. Compare the median and mean expression values using the numpy functions `median` and `mean`.

What can you infer from the difference between these statistics?

4.  In order to work with a more continuous range of expression values, you will need to log-transform the expression values. Because there are zeros in the expression data, you will need to add a pseudo-count (a small value to each expression) to avoid taking the log of zero. Use the numpy function `log2` to transform the data.

Now how do the median and mean transformed expression values compare?

5.  Find the minimum value and maximum value for each gene. You can use the numpy functions `amin` and `amax` to do this. You will need to specify an axis to calculate this statistic across. When you are done, you should have one min value and one max value for each gene. Given that numpy reads axis 0 as rows and axis 1 as columns, think about which axis you need to specify with the `axis=` argument (hint: because the statistic is calculated across the axis specified, the resulting array is the same shape as the input data except excluding the axis specified).

6.  Save the minimum and maximum values for the first 5000 genes to a text file. Load these into R and plot the minimum values (x-axis) by the maximum values (y-axis) in a scatter plot. Make sure to label you axes. It is also helpful to set the point `alpha` to 0.2.

What trends can you see from this plot? How many clusters do you see and can you infer anything about the nature of the expression patterns in each cluster?

7.  Create a copy of the expression array sorted across tissues by expression values (use `numpy.sort`). For each gene, find the difference between the highest expression value and the second highest expression value (referred to below as the `diff_array`). You want to be able to identify genes that are highly tissue-specific.

8.  Create an array of the same shape as the expression array but filled with zeros (you can either use `numpy.zeros_like` or `numpy.zeros` and `data_array.shape`). For each tissue, mark with a one which gene has the highest expression. To do this, you will use the function `numpy.argmax`, which returns the index of the highest value. In order to do this across the whole array instead of using a for loop (for loops are always slow in python, whereas numpy is extremely fast), you can use two list of indices, the first for the row and the second for the column:

```python
high_tissue = numpy.zeros_like(data)
high_tissue[numpy.arange(data.shape[0]), numpy.argmax(data, axis=1)] = 1
```

9.  Set genes with diff_array values < 10 to zero. To do this, you first will create a boolean array from the diff_array using the following syntax (`(diff_array < 10)`). You can then reshape this array to add another axis since it is 1D and the `high_tissue` array is 2D. To reshape it, use the syntax (`(diff_array < 10).reshape(-1, 1)`). Finally, multiply the `high_tissue` array by this reshaped array. The reshaped array will be expanded to match the dimensions of the `high_tissue` array. The result will be a 2D array with ones in positions where highly tissue-specific genes are the highest expressed for a given tissue.

10.  Find which gene-tissue pairs are marked as highly expressed and tissue specific using `numpy.where`. Remember that this will return a tuple (like a list) with two 1D arrays, the first for axis 0 (rows/genes) and the second for axis 1 (columns/tissues). Print out and save each gene-tissue pair. Remember that you have the gene names and tissue names saved and these can be indexed using the appropriate 1D array resulting from the `numpy.where` command.

```python
# E.G. if marker_genes is the 2D array resulting from step 9
high_GTs = numpy.where(marker_genes)
high_genes = high_GTs[0]
high_tissues = high_GTs[1]
# The first gene would be
gene_IDs[high_genes[0]]
# The first tissue would be
tissue_names[high_tissues[0]]
```

Do these genes make sense for the corresponding tissue? You can look up information about each gene at [GeneCards](https://www.genecards.org/).