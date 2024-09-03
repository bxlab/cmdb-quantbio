## Using dictionaries to pull specific samples from GTEx data

For this assignment, you will be looking at how much tissue expression varies across individuals for each of the highly expressed and tissue specific genes. To do this, you will be using the gene-tissue pair results from the morning advanced excercise (which is provided) to pull individual sample expression values for each of these gene-tissue pairs. You will also need the complete GTEx expression data file `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct` and the sample attribute file `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`.

For a description of GTEx, see the [GTEx Portal](https://gtexportal.org/home/aboutAdultGtex).

Because the expression data file is so large, you should create a smaller file to test and debug you script on. You can do this using the command:

```bash
head -n 500 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct > test_data.gct
```

Please submit your answers as a `.py` script with comments explaining the logic of each code block, a `.R` script for plotting the results, and a `.pdf` plot (step #7). (Remember that comment lines start with `#` and are ignored by Python and R). For questions regarding interpretation or discussion of your results, please include your answers as comments interspersed with your code.

1. The first step is to load the gene-tissue pairs from the `gene_tissue.tsv` file. This file contains three columns, geneID, gene name, and tissue, separated by tabs. These data are genes that are highly expressed in a single tissue and the corresponding tissue. Because you will be checking this information later by the geneID, it makes the best sense to create a dictionary keyed by the geneID with the tissue as the value. 


2. Next you will need to figure out which tissue corresponds to which sampleIDs. If you look at the complete gene expression file `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct` (you could use `less` to look at the file), you will see that each column is labeled by a sampleID. Looking at the metadata file `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`, you will see that each sampleID corresponds to a particular subject and tissue. The columns that you are interested are `SAMPID` (the sample ID) and `SMTSD` (the specific tissue label). In the following steps you will be checking things by the tissue name so you will want to create a dictionary using the tissue as the key.
- For each tissue, you will have multiple sample IDs so you will need the dictionary value to be a list that you append sample IDs to. (so you don't need to check if a tissue has already been added to the dictionary, you can use `setdefault`, e.g. `sample_dict.setdefault(tissue, [])`)
- Make sure to look at the structure of the file to know if there are lines that need skipping. In order to skip a line, you need to first open the file (`fs = open(fname)`) and then you can use `readline` to get one line at a time. If you wanted to skip 4 lines, you would call `fs.readline()` 4 times and not save the returned line. You can still read lines in using a `for` loop after using `readline` using the syntax `for line in fs:`. The python will continue reading in lines where you left off with `readline`.


3. You now need to get the list of sampleIDs that are present in the gene expression file `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct`. To do this, you will need the column names (sampleIDs).
- Use `less` to take a look at the first few lines of the file.
- Load the sample IDs (and only the sampleIDs) that correspond to the columns in the GTEx expression data file into a list. Pay attention to which row they are in and which column the sampleIDs start.
- You do not need to load the complete file to get this information, only up to the line containing the column header.


4. Because each sampleID corresponds to a specific tissue and you are only interested in specific tissues for specific genes, it will be useful to know which columns (i.e. which sampleIDs) are relevant for each tissue. Unfortunately, not every sampleID is in the expression data file. To get the relevant columns, the first step is to create a dictionary keyed by each sampleID (from step 3) with the column index as its value.


5. Now you can create the your dictionary of each tissue and which columns in the gene expression file correspond to it. To do this, you will need to create a new dictionary keyed by tissue names with the value being a list of column indexes. To do this, use a nested `for` loop to step through the tissue names followed by stepping through the sampleIDs associated with that tissue (the dictionary you created in step 2). For each sampleID, check if it is present in your sampleID-column index dictionarey from step 4. If it is, add that column index to the correct tissue list in your new dicionary.
- Remember that you can check if an ID is in the sample ID list using a command like `if sample in sample_list:`.

For each tissue type, see how many samples have expression data. Which tissue types have the largest number of samples? The fewest?


6. Now that you know which columns you need for any given tissue, you can load the expression file, keeping only genes that appear in the gene-tissue pair file and only epxression values from sampleIDs that correspond to the tissue of interest for that gene.
- To check if the gene is from the gene-tissue pair file, you can use the keyword `in` to see if a value is in the step 1 dictionary's keys just like a list.
- If the gene is in the gene-tissue pair set, determine which tissue that gene is associated with, get the column indices for that tissue, and pull out only the expression values for the corresponding tissue.
- Unlike lists, numpy arrays can use a list of indices to pull out multiple values at the same time and much faster. With this in mind, when you find a target gene you should convert its expression values into a numpy array and then you can use the index list that you made step 5 to extract the tissue-specific expression values in one step.
- You will potentially have different numbers of expression values for each gene, so saving your data in a numpy array doesn't make sense. Instead, you can use a couple of different approaches. A list with the geneID, tissue, and expressions is one option. A dictionary keyed by the geneID and tissue name with the expression value array as the coresponding value is another.


7. Now that you have the relevant expression values, you need to save the results in a tab-separated file with one line per expression value, its corresponding geneID, and tissue. Your output should have 3 columns. You can do this using a nested `for` loop. The outer `for` loop looks at each gene while the inned `for` loop reads each expression value for that gene.


8. Finally, you can visualize how variable each gene's expression is. Load the data into R and create a violin plot of expression levels broken down by gene (ggplot2's `geom_violin()`).
- For categories, create a combination of tissue names and gene IDs (`dplyr::mutate(Tissue_Gene=paste0(Tissue, " ", GeneID))`)
- You will need to log-transform your data with a psuedo-count of one (you can use `dplyr::mutate` for this step as well)
- Switch the axes for the violin plot so the categories are on the y-axis (`coord_flip()`)
- Make sure to label your axes

Given the tissue specificity and high expression level of these genes, are you surprised by the results? What tissue-specific differences do you see in expression variability? Speculate on why certain tissues show low variability while others show much higher expression variability.


### Advanced exercises

9. For tissues with more than 1 gene, find the correlation between each pairwise combination of gene sample expression levels for that tissue and report the tissue name, number of associated genes, and median correlation (you can use the numpy `corrcoef` to find the Pearson correlation, and because this function returns a correlation matrix, you will need to select item `[0, 1]`).


10. For any tissue with a median correlation greater than 0.4, find and print the gene names and their associated tissue.

Look up the functions of the genes using [GeneCards](https://www.genecards.org/) (you need to remove the period and number following it for gene IDs to look up the genes). Do genes associated with a given tissue appear to be related in function? If so, what function?
