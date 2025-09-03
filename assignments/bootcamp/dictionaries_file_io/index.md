# Assigment 2 - Python dictionaries and file I/O

### Overview
------------

The GTEx gene expression dataset consists of samples taken from many human subjects across a variety of tissues. From this dataset we have identified a set of highly-expressed tissue-specific genes. In this exercise, you will be pulling the expressing data for those genes but only in their primary tissue of expression and investigating how much variability exists across the sample population.

Biological Learning Objectives

- Use metadata to associate data with subject or sample traits
- Use descriptive statistics to make inferences about gene expression heterogeneity

Computational Learning Objectives

- Create and modify a dicitionary
- Use keys to retrieve values stored in a dictionary
- Step through each key and value in a dictionary inside a `for` loop
- Open and read a text file
- Write to a text file

### Instructions
----------------

- Save a copy of this notebook as `~/qbXX-answers/python_dictionaries.ipynb`.
- Fill in answers in the available code/markdown cells below.
- Remember to comment your code to help yourself and us know what each part is intended to do.


```python
# Use this dictionary for identifying genes of interest and their corresponding tissue of interest

gene_tissue = {
    'ENSG00000042832.11': 'Thyroid',
    'ENSG00000091704.9': 'Pancreas',
    'ENSG00000118245.2': 'Testis',
    'ENSG00000122304.10': 'Testis',
    'ENSG00000122852.14': 'Lung',
    'ENSG00000132693.12': 'Liver',
    'ENSG00000134812.7': 'Stomach',
    'ENSG00000135346.8': 'Pituitary',
    'ENSG00000137392.9': 'Pancreas',
    'ENSG00000142515.14': 'Prostate',
    'ENSG00000142615.7': 'Pancreas',
    'ENSG00000142789.19': 'Pancreas',
    'ENSG00000158874.11': 'Liver',
    'ENSG00000162438.11': 'Pancreas',
    'ENSG00000164816.7': 'Small Intestine - Terminal Ileum',
    'ENSG00000164822.4': 'Small Intestine - Terminal Ileum',
    'ENSG00000168925.10': 'Pancreas',
    'ENSG00000168928.12': 'Pancreas',
    'ENSG00000170890.13': 'Pancreas',
    'ENSG00000171195.10': 'Minor Salivary Gland',
    'ENSG00000172179.11': 'Pituitary',
    'ENSG00000175535.6': 'Pancreas',
    'ENSG00000175646.3': 'Testis',
    'ENSG00000182333.14': 'Stomach',
    'ENSG00000187021.14': 'Pancreas',
    'ENSG00000203784.2': 'Testis',
    'ENSG00000204983.13': 'Pancreas',
    'ENSG00000219073.7': 'Pancreas',
    'ENSG00000229859.9': 'Stomach',
    'ENSG00000240338.5': 'Pancreas',
    'ENSG00000254647.6': 'Pancreas',
    'ENSG00000256713.7': 'Stomach',
    'ENSG00000259384.6': 'Pituitary',
    }

```

1. Load sample metadata from "~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" and create a dictionary of sample metadata associating each sample ID with the tissue it was sampled from

    - This file is tab-separated but some values have spaces in their text so make sure to specify `\t` in the `.split()` method
    - There are several columns in the metadata file, but you only need to be concerned with `SAMPID` and `SMTSD`, the sample ID and sample tissue, respectively
    - Your dictionary should use the `SAMPID` as the key and `SMTSD` as the value

```python
metadata_fname = "/Users/cmdb/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

# Create a dictionary to hold the metadata

# Step through each line of the metadata file

    # Remove the newline character from the line you just read in and split it by the tab character

    # Add the data to your metadata dictionary, using the "SAMPID" column value as the dictionary key and the "SMTSD" column value as the dictionary value

```

2. Load and process the headeer from the GTEx gene expression file "~/Data/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"

    - You have code that will open the file and read in the header line
    - Create a list of tissues corresponding to the column names in the header, one tissue per column name (don't worry about which the column is a sample ID or not since column names that aren't sample IDs will just get a "missing" label)
    - For each value in the header, see if it appears in the dictionary of sample ID/tissues that you created in step 1
    - If the value is in your dictionary, add the corresponding tissue name (the value in dictionary for that sample ID key) to your tissue list
    - If the value is not in your dictionary, add "missing" to the tissue list

```python
# Because this file is a gzip-compressed file, you will be using the `gzip` library that comes with python
import gzip

expression_fname = "/Users/cmdb/Data/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"

# Open the file using the gzip library
fs = gzip.open(expression_fname)

# The file starts with two lines of information before the header, so these are skipped using 2 calls of `.readline()` without keeping the returned data
_ = fs.readline()
_ = fs.readline()

# Gzipped text files are read in as byte strings, not regular strings, so you will see that `.decode()` is included after `.readline()` to convert the input into a string
header = fs.readline().decode().rstrip().split("\t")

# Make sure to close the file when you are done reading it
fs.close()

# Create a list to hold the tissue names
tissues = []

# Step through each value in the header

    # Check if the value has an entry in the metadata dictionary, adding either the tissue corresponding to that entry or "missing" to your tissues list

```

3. Create a list of expression values for each gene in the `gene_tissue` dictionary (given at the start of the exercise), keeping only expression values from the tissue of interest for that gene

    - Since the file is already open, you can keep reading lines from it in, picking up from after the last line read
    - You will be log-transforming the expression data using the `log2` function from the built-in library `math`
    - To deal with zeros in expression data, add one before log-transforming expression data

```python
# Import the log2 function from the math library to transform the expression counts
from math import log2

# Create a dictionary to hold each gene's expression data
gene_expression = {}

# Because this is a large file, for testing purposes you can use only a small portion of the file to test your code
line_counter = 0

# Step through each line of the file
for line in gzip.open(expression_fname):

    ############# REMOVE THIS AND RERUN ONCE YOUR CODE IS WORKING ####################
    # For debugging, use only the first 500 lines of the file.
    line_counter += 1
    if line_counter == 500:
        break
    ############# REMOVE THIS AND RERUN ONCE YOUR CODE IS WORKING ####################
    

    # Convert the line into a string from a byte string, remove the newline character at the end, and split it by tab characters
    fields = line.decode().rstrip().split("\t")

    # Check if this line contains data for a gene we are interested in. If not, skip it

    # Determine which tissue we care about for this gene

    # Create a list in the gene expression dictionary to hold the gene's data
    gene_expression[fields[0]] = []

    # Look at each column position, one by one
    for i in range(len(tissues)):


        # For each column, see if it is from the tissue we care about (using the tissues list). If so, add that field to our gene expression list, first transforming it with log2(1 + float(field[i])) 

```

4. For each gene, calculate the mean, standard deviation, and relative standard deviation and write them to a results file

    - You will be using the `mean` and `stdev` functions from the built-in `statistics` library. Both functions accept a list as their input
    - Format your results for each gene into a string and write that string to your results file
    - If you want to organize your results better (completely optional!!), you can save your results in a list and sorting them by tissue type before writing them to a file

```python
# Load the functions to caculcate the mean and standard deviation from a list
from statistics import stdev, mean

# Open a file to write the results to

# Step through each gene and its expression value list

    # Calculate the mean and standard deviation of expression

    # Calculate the relative standard deviation (coefficient of variation, std / mean)

    # Write the tissue namne, gene name, mean, standard deviation, and relative standard deviation to your results file

# Close your results file

```

5. Looking at your results, answer the following questions

    - Which tissues have low relative variability (rv < 0.2), medium relative variability (0.2 < rv < 0.4), or high relative variability (rv > 0.4)?
    - Do each of the genes for a given tissue show the same level of relative variability?
    - Propose an explanation for differences in relative variability you see between tissues based on what you know about their functions.
    - Were any of the results surpising to you?
