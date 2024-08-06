## Using dictionaries to pull specific samples from GTEx data

For a description of GTEx, see the [GTEx Portal](https://gtexportal.org/home/aboutAdultGtex).

You will be using the gene-tissue pair results from the morning excercise to pull individual sample expression values for each of these gene-tissue pairs. To do this you will also need the complete GTEx expression data file `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct` and the sample attribute file `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`.

Please submit your answers as a `.py` script with comments explaining the logic of each code block, a `.R` script for plotting the results, and a `.pdf` plot (step #7). (Remember that comment lines start with `#` and are ignored by Python and R). For questions regarding interpretation or discussion of your results, please include your answers as comments interspersed with your code.

1.  Load the geneID and tissue pairs into a dictionary such that the geneID is the key and the tissue is the value.

```python
# Load gene-tissue pair data
GT_dict = {}
for line in open(gene_tissue_fname):
    gene_ID, gene_name, tissue = line.rstrip().split("\t")
    GT_dict[gene_ID] = tissue
```

2.  Load the sample attributes into a dictionary. Make sure to look at the structure of the file to know if there are lines that need skipping. The columns that you are interested are `SAMPID` (the sample ID) and `SMTSD` (the specific tissue label). You will want to use the tissue as the key. For each tissue, you will have multiple sample IDs so you will need the value to be a list that you append sample IDs to.

```python
# Load sample IDs for tissue association
sample_dict = {}
fs = open(attribute_fname)
_ = fs.readline()
for line in fs:
    fields = line.rstrip().split("\t")
    sample = fields[0]
    tissue = fields[6]
    sample_dict.setdefault(tissue, [])
    sample_dict[tissue].append(sample)
```

3.  Load the sample IDs that correspond to the columns in the GTEx expression data file. Pay attention to which row they are in and where the sample IDs start.

```python
# Load expression data
genes = []
data = []
fs = open(data_fname)
# Skip header lines
_ = fs.readline()
_ = fs.readline()
# Save sample IDs for colum  identification
sample_IDs = fs.readline().rstrip().split("\t")[2:]
```

4.  Not every sample ID in the attribute file is in the expression data file so you will need to convert the list of sample IDs in the dictionary you created in step 2 into a list of column indices. You can create a new dictionary keyed by tissue type with the index list as the values or just replace the list of sample IDs with the list of indices in your current sample dictionary.
- Remember that you can check if an ID is in the sample ID list using a command like `if sample in sample_list:`
- You can get the position of a specific value in a list using the list method `index`, e.g. `position = sample_list.index(sample)`.

```python
# Convert sample IDs to column indices in sample dict
for tissue in sample_dict.keys():
    indices = []
    for sample_ID in sample_dict[tissue]:
        # Check that sample_ID is in data sample_IDs
        if sample_ID in sample_IDs:
            indices.append(sample_IDs.index(sample_ID))
    sample_dict[tissue] = numpy.array(indices)
    print(tissue, len(indices))
```

For each tissue type, see how many samples have expression data. Which tissue types have the largest number of samples? The fewest?

**Blood, skin, and other externally-accessible tissues have the most samples which brain, internal organs, and particularly female reproductive organ tissues have many fewer samples.**

5.  Now that you know which columns you need for any given tissue, load the expression file, keeping only genes that appear in the gene-tissue pair file. You can use the keyword `in` to see if a value is in a dictionary's keys just like list. If the gene is in the gene-tissue pair set, pull out only the expression values for the corresponding tissue.
- Unlike lists, numpy arrays can use a list of indices to pull out multiple values at the same time and much faster. With this in mind, when you find a target gene you should convert its expression values into a numpy array and then you can use the indice list that you made step #4 to extract the tissue-specific expression values.

```python
# Load data for each gene
for line in fs:
    fields = line.rstrip().split('\t')
    gene = fields[0]
    # If not a marker gene, skip
    if gene not in GT_dict:
        continue
    # Record gene to keep track of order
    genes.append(gene)
    # Convert expressions to a 1D array
    expr = numpy.array(fields[2:], dtype=float)
    # Get tissue associated with marker gene
    tissue = GT_dict[gene]
    # Get tissue sample indices
    sample_indices = sample_dict[tissue]
    # Keep tissue expressions
    data.append(expr[sample_indices])
fs.close()
```
6.  Save the results in a tab-separated file with one line per expression value, its corresponding geneID, and tissue.

```python
# Save expressions by tissue
output = open(out_fname, 'w')
output.write(f"GeneID\tTissue\tExpr\n")
for i in range(len(genes)):
    gene = genes[i]
    tissue = GT_dict[gene]
    for j in range(len(data[i])):
        output.write(f"{gene}\t{tissue}\t{data[i][j]}\n")
output.close()
```

7.  Load the data into R and create a violin plot of expression levels broken down by gene.
- You will need to log-transform your data with a psuedo-count of one
- For categories, create a combination of tissue names and gene IDs
- Switch the axes for the violin plot so the categories are on the y-axis
- Make sure to label your axes

```R
library(readr)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

data = readr::read_delim(args[1])

# Add log-transformed data column
data = data %>%
    dplyr::mutate(Log2_Expr=log2(Expr + 1))

# Add tissue-gene data column
data = data %>%
    dplyr::mutate(Tissue_Gene=paste0(Tissue, " ", GeneID))

# Create plot
p = ggplot(data, aes(x=Tissue_Gene, y=Log2_Expr)) +
    geom_violin() +
    coord_flip() +
    xlab("Tissue + Gene") +
    ylab("Log2 Expression")

# Save plot
pdf(args[2])
print(p)
dev.off()
```

Given the tissue specificity and high expression level of these genes, are you surprised by the results? What tissue-specific differences do you see in expression variability? Speculate on why certain tissues show low variability while others show much higher expression variability.

**There is a wide range of variability in gene expression between individuals but the amount of variability is quite consistent across genes found in the same tissue. The testis, thyroid, and pancreas show much lower levels of variability compared to the stomach and small intestine, which makes sense given the crucial roles of the former. What is surprising is the variability seen in the pituitary-associated genes.**