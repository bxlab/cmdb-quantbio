## Using Numpy matrices and indexing to find marker genes in GTEx data

For a description of GTEx, see the [GTEx Portal](https://gtexportal.org/home/aboutAdultGtex).

Certain genes have high tissue specificity while some tissues are highly enriched in certain proteins. In this assignment, you will be find the intersection between these two, identifying genes with a ten-fold enrichment in a single tissue and that are the highest expressed gene in their enriched tissue. To do this, you will be using a summary dataset that has the median expression levels across samples for each gene and tissue, named `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct`.

Please submit your answers as a `.py` script with comments explaining the logic of each code block and for questions regarding interpretation or discussion of your results, please include your answers as comments interspersed with your code or a separate text file (Remember that comment lines start with `#` and are ignored by Python and R). For advanced exercises, please also include the `.tsv` file with gene-tissue pairs. 


#### 1. The first step is to load the expression data into a nested list, the tissue names into a list, the gene IDs into a list, and the gene names into a list.
- These are all in the same file so you will need to parse them as you read the file.
- Load the GTEx file using the standard file `open` function discussed in previous sections.
- In order to load the data, first look at the file formatting using `less GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct` on the command line. Note that there are lines at the beginning that you want to skip.You can use the filestream method `readline()` to read in one line at a time. This will allow you to read in header lines without saving the results. If you wanted to skip 4 lines, you would call `fs.readline()` 4 times and not save the returned line. You can still read lines in using a `for` loop after using `readline` using the syntax `for line in fs:`. Python will continue reading in lines from where you left off with `readline`.
- There are four types of information you will want to capture: tissue names, gene IDs, gene names, and expression levels. Each field is separated by a `\t` (a tab), so make sure to specify that when using the `split()` method.


#### 2.  Because nested lists are clunky and difficult to work with, you will  be converting the data into numpy arrays.
- Convert each of your data lists into a numpy array.
- In the case of text names, numpy will determine what datatype is appropriate. For your expression data, you will need to specify the data type as `float` using the argument `dtype=float`.

   ***Why do you need to tell numpy what type of data the expression data is but not for the other data lists?***


#### ~~3. Let's check the mean expression values for the first 10 genes using an approach you have already learned, nested `for` loops.~~
~~For the first 10 genes (rows, axis 0), use a nested `for` loop to find the mean gene expression and print them. Think about which axis needs to be the outer `for` loop and which the inner.~~


#### 4. Now let's see how numpy arrays make working with data more streamlined by calculating the same mean expression values for the first 10 genes but using the build-in numpy function `mean` and printing them.
- You will need to use the `axis` argument and set it equal to the axis corresponding to tissues since you want to find the average across tissues.

    ***Do they match the means you found with the nested `for` loop approach?***


#### 5. To get a sense of the spread of data, calculate and compare the median and mean expression value of the entire dataset.
- so that the statistics are calculated across all values, do **not** use the `axis` argument.
- The numpy functions for these statistics are `median` and `mean`.

    ***What can you infer from the difference between these statistics?***


#### 6. In order to work with a more normalized range of expression values, let's apply a log-transformation to the data and check the mean and median values again.
- In order to keep the values interpretable, use a log2 (base 2) transformation (the function `numpy.log2`).
- Because there are zeros in the expression data, you will need to add a pseudo-count (a small value to each expression) to avoid taking the log of zero. In this case, use a pseudo-count of one.

    ***Now how do the median and mean transformed expression values compare to each other? To the non-transformed values?***


#### 7. Now let's find the expression gap for each gene between their highest and next highest expression level to identify highly tissue specific genes.
- To do this, you will need to create a copy of the expression array since you will be altering the array to find this information. In order to make a copy of an array, you can use the function `numpy.copy`.
- Sort this copy across tissues by expression values (use `numpy.sort`). You will need to specify an axis to sort along. Given that numpy reads axis 0 as rows and axis 1 as columns, think about which axis you need to specify with the `axis=` argument.
- For each gene, find the difference between the highest expression value and the second highest expression value (referred to below as the `diff_array`). You want to be able to identify genes that are highly tissue-specific.
- Hint: remember that you can use negative number indices to work from the last item in a list/array backward


#### 8. Finally, using the expression difference array you just created, you can now identify genes that show high single-tissue specificity as defined as a difference of at least 10 (~1000-fold difference since the data are log2-transformed).
- Print the number of genes whose difference between the highest and second highest tissue expression is greater than 10.


### Advanced exercises

#### 9. You now need to figure out for each tissue which gene is the most highly expressed.
- To do this, create an array of the same shape as the expression array but filled with zeros (you can either use `numpy.zeros_like` or `numpy.zeros` and `data_array.shape`).
- For each tissue, mark with a one which gene has the highest expression. To do this, you will use the function `numpy.argmax`, which returns the index of the highest value.
- In order to do this across the whole array instead of using a `for` loop (`for` loops are always slow in python, whereas numpy is extremely fast), you can use two list of indices, the first for the row and the second for the column:

```python
row_index = numpy.arange(data.shape[0])
col_index = numpy.argmax(...fill in yourself...)
```

#### 10. Next, you will want to eliminate genes that do not have a high tissue specificity as determined in step 8 by setting genes with diff_array values <= 10 to zero.
- To do this, you first will create a boolean array from `diff_array` using the following syntax (`(diff_array < 10)`).
- You can then reshape this array to add another axis since it is 1D and the `high_tissue` array is 2D. To reshape it, use the syntax (`(diff_array < 10).reshape(-1, 1)`).
- Finally, multiply the `high_tissue` array by this reshaped array. The reshaped array will be expanded to match the dimensions of the `high_tissue` array. The result will be a 2D array with ones in positions where highly tissue-specific genes are the highest expressed for a given tissue.


#### 11. You should now have an array of the same shape as the original expression array, but with ones for gene-tissue pairs with high expression and specificity and zeros everywhere else.
- Find which gene-tissue pairs are marked as highly expressed and tissue specific (i.e. have a value of one) using `numpy.where`.
- Remember that this will return a tuple (like a list) with two 1D arrays, the first for axis 0 (rows/genes) and the second for axis 1 (columns/tissues).
- Print out and save each gene-tissue pair. Remember that you have the gene names and tissue names saved and these can be indexed using the appropriate 1D array resulting from the `numpy.where` command.

    ***Do these genes make sense for the corresponding tissue? You can look up information about each gene at [GeneCards](https://www.genecards.org/) (you need to remove the period and number following it for gene IDs to look up the genes as the number after the last period is a version number).***
