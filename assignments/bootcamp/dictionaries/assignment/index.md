## Using dictionaries to pull specific samples from GTEx data

For a description of GTEx, see the [GTEx Portal](https://gtexportal.org/home/aboutAdultGtex).

You will be using the gene-tissue pair results from the morning advanced excercise to pull individual sample expression values for each of these gene-tissue pairs. To do this you will also need the complete GTEx expression data file `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct` and the sample attribute file `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`.

Because the expression data file is so large, you should create a smaller file to test and debug you script on. You can do this using the command `head -n 500 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct > test_data.gct`.

Please submit your answers as a `.py` script with comments explaining the logic of each code block, a `.R` script for plotting the results, and a `.pdf` plot (step #7). (Remember that comment lines start with `#` and are ignored by Python and R). For questions regarding interpretation or discussion of your results, please include your answers as comments interspersed with your code.

1.  Using the `gene_tissue.tsv` file, load the geneID and tissue pairs into a dictionary such that the geneID is the key and the tissue is the value. This file contains three columns, geneID, gene name, and tissue, separated by tabs. These data are genes that are highly expressed in a single tissue and the corresponding tissue.


2.  Using the file `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`, load the tissue names and the corresponding sample IDs into a dictionary. Make sure to look at the structure of the file (you could use `less` to look at the file) to know if there are lines that need skipping. The columns that you are interested are `SAMPID` (the sample ID) and `SMTSD` (the specific tissue label). You will want to use the tissue as the key. For each tissue, you will have multiple sample IDs so you will need the dictionary value to be a list that you append sample IDs to. (so you don't need to check if a tissue has already been added to the dictionary, you can use `setdefault`, e.g. `sample_dict.setdefault(tissue, [])`)
- In order to skip a line, you need to first open the file (`fs = open(fname)`) and then you can use `readline` to get one line at a time. If you wanted to skip 4 lines, you would call `fs.readline()` 4 times and not save the returned line. You can still read lines in using a `for` loop after using `readline` using the syntax `for line in fs:`. The python will continue reading in lines where you left off with `readline`.


3.  In order to know which samples are in the expression file `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct`, you will need the column names. Use `less` to take a look at the first few lines of the file. Load the sample IDs (and only the sample IDs) that correspond to the columns in the GTEx expression data file into a list. Pay attention to which row they are in and which column the sample IDs start.


4.  Not every sample ID in the attribute file is in the expression data file so you will need to filter and convert the list of sample IDs in the dictionary you created in step 2 into a list of column indices. You can create a new dictionary keyed by tissue type with the index list as the value or just replace the list of sample IDs with the list of indices in your current sample dictionary. This will require a nested `for` loop. The outside `for` loop should step through each tissue in your sample dictionary. The inner `for` loop should look at each sample ID in the sample dictionary value corresponding to the tissue (the current value of the outer `for` loop). 
- Remember that you can check if an ID is in the sample ID list using a command like `if sample in sample_list:`
- You can get the position of a specific value in a list using the list method `index`, e.g. `position = sample_list.index(sample)`.

For each tissue type, see how many samples have expression data. Which tissue types have the largest number of samples? The fewest?


5.  Now that you know which columns you need for any given tissue, load the expression file, keeping only genes that appear in the gene-tissue pair file. You can use the keyword `in` to see if a value is in a dictionary's keys just like a list. If the gene is in the gene-tissue pair set, determine which tissue that gene is associated with, get the column indices for that tissue, and pull out only the expression values for the corresponding tissue.
- Unlike lists, numpy arrays can use a list of indices to pull out multiple values at the same time and much faster. With this in mind, when you find a target gene you should convert its expression values into a numpy array and then you can use the indice list that you made step #4 to extract the tissue-specific expression values in one step.


6.  Save the results in a tab-separated file with one line per expression value, its corresponding geneID, and tissue. Your output should have 3 columns. You can do this using a nested `for` loop. The outer `for` loop looks at each gene while the inned `for` loop reads each expression value for that gene.


7.  Load the data into R and create a violin plot of expression levels broken down by gene (ggplot2's `geom_violin()`).
- For categories, create a combination of tissue names and gene IDs (`dplyr::mutate(Tissue_Gene=paste0(Tissue, " ", GeneID))`)
- You will need to log-transform your data with a psuedo-count of one (you can use `dplyr::mutate` for this step as well)
- Switch the axes for the violin plot so the categories are on the y-axis (`coord_flip()`)
- Make sure to label your axes

Given the tissue specificity and high expression level of these genes, are you surprised by the results? What tissue-specific differences do you see in expression variability? Speculate on why certain tissues show low variability while others show much higher expression variability.


### Advanced exercises

8. For tissues with more than 1 gene, find the correlation between each pairwise combination of gene sample expression levels for that tissue and report the tissue name, number of associated genes, and median correlation (you can use the numpy `corrcoef` to find the Pearson correlation, and because this function returns a correlation matrix, you will need to select item `[0, 1]`).


9. For any tissue with a median correlation greater than 0.4, find and print the gene names and their associated tissue.

Look up the functions of the genes using [GeneCards](https://www.genecards.org/) (you need to remove the period and number following it for gene IDs to look up the genes). Do genes associated with a given tissue appear to be related in function? If so, what function?
