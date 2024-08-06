## Using dictionaries to pull specific samples from GTEx data

For a description of GTEx, see the [GTEx Portal](https://gtexportal.org/home/aboutAdultGtex).

You will be using the gene-tissue pair results from the morning excercise to pull individual sample expression values for each of these gene-tissue pairs. To do this you will also need the complete GTEx expression data file `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct` and the sample attribute file `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`.

Please submit your answers as a `.py` script with comments explaining the logic of each code block, a `.R` script for plotting the results, and a `.pdf` plot (step #7). (Remember that comment lines start with `#` and are ignored by Python and R). For questions regarding interpretation or discussion of your results, please include your answers as comments interspersed with your code.

1.  Load the geneID and tissue pairs into a dictionary such that the geneID is the key and the tissue is the value.

2.  Load the sample attributes into a dictionary. Make sure to look at the structure of the file to know if there are lines that need skipping. The columns that you are interested are `SAMPID` (the sample ID) and `SMTSD` (the specific tissue label). You will want to use the tissue as the key. For each tissue, you will have multiple sample IDs so you will need the value to be a list that you append sample IDs to.

3.  Load the sample IDs that correspond to the columns in the GTEx expression data file. Pay attention to which row they are in and where the sample IDs start.

4.  Not every sample ID in the attribute file is in the expression data file so you will need to convert the list of sample IDs in the dictionary you created in step 2 into a list of column indices. You can create a new dictionary keyed by tissue type with the index list as the values or just replace the list of sample IDs with the list of indices in your current sample dictionary.
- Remember that you can check if an ID is in the sample ID list using a command like `if sample in sample_list:`
- You can get the position of a specific value in a list using the list method `index`, e.g. `position = sample_list.index(sample)`.

For each tissue type, see how many samples have expression data. Which tissue types have the largest number of samples? The fewest?

5.  Now that you know which columns you need for any given tissue, load the expression file, keeping only genes that appear in the gene-tissue pair file. You can use the keyword `in` to see if a value is in a dictionary's keys just like list. If the gene is in the gene-tissue pair set, pull out only the expression values for the corresponding tissue.
- Unlike lists, numpy arrays can use a list of indices to pull out multiple values at the same time and much faster. With this in mind, when you find a target gene you should convert its expression values into a numpy array and then you can use the indice list that you made step #4 to extract the tissue-specific expression values.

6.  Save the results in a tab-separated file with one line per expression value, its corresponding geneID, and tissue.

7.  Load the data into R and create a violin plot of expression levels broken down by gene.
- You will need to log-transform your data with a psuedo-count of one
- For categories, create a combination of tissue names and gene IDs
- Switch the axes for the violin plot so the categories are on the y-axis
- Make sure to label your axes

Given the tissue specificity and high expression level of these genes, are you surprised by the results? What tissue-specific differences do you see in expression variability? Speculate on why certain tissues show low variability while others show much higher expression variability.