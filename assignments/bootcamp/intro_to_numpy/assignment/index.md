## Using Numpy matrices and indexing to find marker genes in GTEx data

For a description of GTEx, see the [GTEx Portal](https://gtexportal.org/home/aboutAdultGtex).

You will be using a summary dataset that has the median expression levels across samples for each gene and tissue, named `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct`.

Please submit your answers as a `.py` script with comments explaining the logic of each code block and for questions regarding interpretation or discussion of your results, please include your answers as comments interspersed with your code or a separate text file (Remember that comment lines start with `#` and are ignored by Python and R). For advanced exercises, please also include the `.tsv` file with gene-tissue pairs. 


1. Load the GTEx file using the standard file `open` function discussed in previous sections.
- In order to load the data, first look at the file formatting using `less GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct` on the command line. Note that there are lines at the beginning that you want to skip.
- You can use the filestream method `readline()` to read in one line at a time. This will allow you to read in header lines without saving the results. If you wanted to skip 4 lines, you would call `fs.readline()` 4 times and not save the returned line. You can still read lines in using a `for` loop after using `readline` using the syntax `for line in fs:`. The python will continue reading in lines where you left off with `readline`.
- There are four types of information you will want to capture: tissue names, gene IDs, gene names, and expression levels. Each field is separated by a `\t`, so make sure to specify that when using the `split()` method.
- You should end up with a list of tissue names, a list of gene IDs, a list of gene names, and a nested list of expression levels


2.  Convert each of your data lists into a numpy array. In the case of text names, numpy will determine what datatype is appropriate. For your expression data, you will need to specify the data type as `float` using the argument `dtype=float`.

Why do you need to tell numpy what type of data the expression data are but not for the other data lists?


3. For the first 10 genes (rows, axis 0), use a nested `for` loop to find the mean gene expression and print them. Think about which axis needs to be the outer `for` loop and which the inner.


4. Now use the build-in numpy function `mean` to find the mean expression for the first 10 genes and print them (you will need to use the `axis` argument and set it equal to the axis corresponding to tissues since you want to fund the average across tissues).

Do they match the means you found with the nested `for` loop approach?


5. Print out the median for the entire set of expression values (i.e. do **not** use the `axis` argument) and the mean across all expression values using the numpy functions `median` and `mean`.

What can you infer from the difference between these statistics?


6.  In order to work with a more continuous range of expression values, you will need to log-transform the expression values. Because there are zeros in the expression data, you will need to add a pseudo-count (a small value to each expression) to avoid taking the log of zero. Use the numpy function `log2` to transform the data.

Now how do the median and mean transformed expression values compare?


7.  Create a copy of the expression array. Sort this copy across tissues by expression values (use `numpy.sort`). You will need to specify an axis to sort along. Given that numpy reads axis 0 as rows and axis 1 as columns, think about which axis you need to specify with the `axis=` argument. For each gene, find the difference between the highest expression value and the second highest expression value (referred to below as the `diff_array`). You want to be able to identify genes that are highly tissue-specific.
- Hint: remember that you can use negative number indices to work from the last item in a list/array backward


8. Find and print the number of genes whose difference between the highest and second highest tissue expression is greater than 10 (~1000-fold difference since the data are log2-transformed).


### Advanced exercises

9.  Create an array of the same shape as the expression array but filled with zeros (you can either use `numpy.zeros_like` or `numpy.zeros` and `data_array.shape`). For each tissue, mark with a one which gene has the highest expression. To do this, you will use the function `numpy.argmax`, which returns the index of the highest value. In order to do this across the whole array instead of using a for loop (for loops are always slow in python, whereas numpy is extremely fast), you can use two list of indices, the first for the row and the second for the column:


10.  Set genes with diff_array values < 10 to zero. To do this, you first will create a boolean array from the diff_array using the following syntax (`(diff_array < 10)`). You can then reshape this array to add another axis since it is 1D and the `high_tissue` array is 2D. To reshape it, use the syntax (`(diff_array < 10).reshape(-1, 1)`). Finally, multiply the `high_tissue` array by this reshaped array. The reshaped array will be expanded to match the dimensions of the `high_tissue` array. The result will be a 2D array with ones in positions where highly tissue-specific genes are the highest expressed for a given tissue.


11.  Find which gene-tissue pairs are marked as highly expressed and tissue specific using `numpy.where`. Remember that this will return a tuple (like a list) with two 1D arrays, the first for axis 0 (rows/genes) and the second for axis 1 (columns/tissues). Print out and save each gene-tissue pair. Remember that you have the gene names and tissue names saved and these can be indexed using the appropriate 1D array resulting from the `numpy.where` command.

Do these genes make sense for the corresponding tissue? You can look up information about each gene at [GeneCards](https://www.genecards.org/).
