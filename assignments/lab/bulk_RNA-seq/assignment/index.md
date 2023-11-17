# Bulk RNA-seq and differential expression analysis

## Assignment Overview

Today, you will be examining gene expression data from the [Genotype-Tissue Expression (GTEx) Project](https://doi.org/10.1126/science.aaz1776)--a prominent resource in the field of human genetics comprising whole-genome and RNA-seq data from >800 post-mortem individuals across >50 tissues. Not all individuals donated all tissues, but there are more than 17,000 total RNA-seq samples, offering detailed insight into gene expression patterns across tissues and individuals, as well as the DNA sequence variants that modulate these differences.

One extremely common use of RNA-seq data is to compare patterns of gene expression across different conditions (i.e., predictor variables), which may include features such as sex, tissue type, genotype (e.g., wild-type versus mutant), drug treatment, time, etc. Such statistical analyses are termed "differential expression" tests, as we effectively loop over each gene and ask whether expression of that gene is correlated with the relevant predictor variable. In some study designs, you may want to "model out" or account for the effects of certain predictor variables (termed "covariates") while focusing attention on other predictors. From the statistical perspective, covariates are typically treated in the same way as other predictor variables, we simply place less focus on their effects in downstream analysis and interpretation as they may be less relevant to our key hypotheses.

Unfortunately, the most popular software packages for performing differential expression analysis (edgeR and DESeq2) were written as libraries for the R programming language, so until this year, it was not possible to introduce these packages in QuantBio Lab. As of this year, DESeq2 has been ported to Python and released as a package called "PyDESeq2". While this provides key functions of DESeq2, the package is still poorly documented, so I would like to start by motivating this exercise with a simpler (but flawed) approach for differential expression analysis based on simple linear regression. Remember that DESeq2/PyDESeq2 are also linear regression(!), but with some additional bells and whistles that I reviewed in the lecture.<br><br>

## Data

The data you'll be using today is derived from the paper mentioned above, focusing on the RNA-seq samples from whole blood (755 total individuals). They were downloaded directly from the [GTEx portal](https://gtexportal.org/home/downloads/adult-gtex#bulk_tissue_expression) and slightly reformatted to save you some time on tedious data wrangling. Download from the Dropbox links below:

[Subject-level metadata](https://www.dropbox.com/scl/fi/zidlbn4rlvyv43k022mmn/gtex_metadata.txt?rlkey=j6aidakljr0739tnnzvpbg0gn&dl=0)

[Gene expression matrix](https://www.dropbox.com/scl/fi/7iengpyrevd356dfq53pg/gtex_whole_blood_counts_formatted.txt?rlkey=l5h12cyher33kkzlrwi4qwf8g&dl=0)<br><br>

## Exercises

There are three exercises in this assignment:

1. Manually perform regression testing for differential expression
2. Perform the same tests using the `pydeseq2` library
3. Data visualization 

Before you do anything else, create a Python script for this assignment. Everything you'll need to do for this assignment will be done in this script and submitted via GitHub.<br><br>

### Exercise 1: Perform a "homemade" test for differential expression between the sexes

You're curious whether any genes are differntially expressed between sexes (in whole blood). Before you do the DE analysis the "right" way with DESeq2, you're first going to perform your own naive "homemade" analysis to get an idea of what DESeq2 is doing behind the scenes at a basic level. This will involve regressing expression quantifications from GTEx onto sample sex, and identifying genes with significant differences in expression between sexes.<br><br>

#### **Step 1.1**: Loading data and importing libraries

Because the data are rectangular data of mixed type, you can use `pandas` to read the data into a data frame. You will also need `statsmodels` to perform your "homemade" differential expression test, followed by PyDESeq2 to apply the more sophisticated test. You'll also need `numpy` to do some data transformations. Here is code to load those up, along with the data. You will need to modify the file paths to direct it to wherever you stored the data.

```
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats import multitest
from pydeseq2 import preprocessing
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# read in data
counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

# read in metadata
metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)
```
<br>

#### **Step 1.2**: Normalization

Before running your analysis, you will first use PyDESeq2 to perform normalization across samples to account for differences in sequencing depth and RNA composition (i.e., some highly expressed genes "eating up" lots of reads and distorting the patterns). This will produce a "normalized counts" matrix that we can use for our own test.

```
counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]
```

In addition to normalizing the counts, you will also want to log them. OLS regression (the type that you will be using for your homemade approach) assumes the data is normally distributed. Unfortunately, gene expression counts definitely are NOT. Logging the counts will make the data closer to a normal distribution, however. We can log the normalized expression counts using `numpy`.

```
counts_df_normed = np.log2(counts_df_normed + 1)
```

<span style="color:red;font-weight:bold">IMPORTANT NOTE:</span> DESeq2 doesn’t actually use these normalized counts as input, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). So when you proceed to Exercise 2,  you will use `counts_df`, not `counts_df_normed`!<br><br>

#### **Step 1.3**: Create a design matrix

To make the regression easier, you're going to put all of your data (`counts_df_normed` and `metadata`) into a single data frame. You can do this using `pd.concat()`

```
full_design_df = pd.concat([counts_df_normed, metadata], axis=1)
```

Take a look at how `full_design_df` looks. You should have one row per sample, with logged expression values 
for all genes, and metadata for that sample.<br><br>

#### **Step 1.4**: Running regression for a single gene

To get started, just run your regression for a single gene. You can start with the first gene in the dataframe: DDX11L1.

You will use statsmodels to perform the statistical test, using the gene ID as the response variable and sex as the predictor variable. Note that the format of some of the gene IDs in your data violate the formatting requirements of `smf.ols`. You can get around this by wrapping the response variable in `Q()` in your regression formula. Note as well that sex here is encoded as 1 and 2, where 1 refers to males and 2 refers to females.

```
model = smf.ols(formula = 'Q("DDX11L1") ~ SEX', data=full_design_df)
results = model.fit()
```

Now, extract and examine the slope and p-value.

```
slope = results.params[1]
pval = results.pvalues[1]
```
<br>

#### **Step 1.5**: Extend this test to all genes

Write a for-loop in Python to extend this test to all genes in your matrix. For each gene that you test, store the slopes and p-values, along with the associated gene names in a useful data structure (pandas DataFrame, arrays, ..., ?).

NOTE: This step will probably take a few minutes to run. We recommend outputing your results to a file, and then reading it back in to your script so that you don't have to re-run this whole process each time you run your python script.

After running the analysis on all the genes, use the Benjamini-Hochberg procedure to determine which genes are sex-differentially expressed at an FDR of 10%. Statsmodels has a great function to do this correction for you: [statsmodels.stats.multitest.fdrcorrection](https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.fdrcorrection.html)

Write the list of genes that are differentially expressed at a 10% FDR to a file. You will be uploading this file with your assignment.<br><br>

### Exercise 2: Repeat differential expression testing with PyDESeq2

Now that you've run your homemade analysis to get a feel for what the analysis is doing, you want to run it the "right" way with DESeq2.

First load the data into a PyDESeq2 object:
```
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors="SEX",
    n_cpus=4,
)
```

Then apply the differential expression test and extract the results:

```
dds.deseq2()
stat_res = DeseqStats(dds)
stat_res.summary()
results = stat_res.results_df
```

Note that the `padj` column of `results` reports the FDR-adjusted p-value (i.e., "q-value"). The rows with a `padj` < 0.1 are the genes that are differentially expressed at an FDR of 10%. A lot of rows will have missing `padj` values. This is because DESeq2 does a filtering step to remove genes that it thinks have a low probability of being DE, and only performs FDR correction for the set of genes it does not filter out. You can ignore genes with missing `padj` values.

Compare the list of genes that ARE differentially expressed at a 10% FDR to those you identified in Step 1.4. What is the percentage of overlap? Compute this percentage in your code as a "Jaccard index", which is defined as the intersection divided by the union: `((number of genes that were significant in steps 1 and 2) / (number of genes that were significant in steps 1 or 2)) * 100%`. Record this in your `README.md` file for this assignment, which you will upload with the assignment.<br><br>

### Exercise 3: Visualization

Use `matplotlib` to create a "Volcano" plot depicting your differential expression results from Exercise 2. A volcano plot is a scatter plot, where the x-axis shows the `log2FoldChange` and the y-axis shows the -log10(pvalue).

Highlight the genes that are significant at a 10% FDR and for which the absolute value of the log2FoldChange is greater than 1 in a separate color.

Output this plot to a `.png` that you will upload with your assignment.<br><br>

## Submission

1. Python script to run DE analysis **(6 points total)**
  * Implementation of manual DE test (Steps 1.1-1.4) **(1 point)**
  * Code to perform DE test on all genes (Step 1.5) **(1 point)**
  * Code for FDR correction (Step 1.5) **(1 point)**
  * Running PyDESeq2 on all genes (Exercise 2) **(1 point)**
  * Code for percent overlap in results between methods (Exercise 2) **(1 point)**
  * Code for volcano plot (Exercise 2) **(1 point)**
2. `README.md` file with answers to questions **(1 point total)**
  *  Percent overlap between manual testing and PyDESeq2 from Exercise 2 **(1 point)**
3. Output text files **(2 points total)**
  *  List of differentially expressed transcripts (10% FDR) from Step 1.5 **(1 point)**
  *  List of differentially expressed transcripts (10% FDR) from Exercise 2 **(1 point)**
4. Plots **(1 point total)**
  * Nicely formatted and labelled volcano plot from Exercise 2**(1 point)**

**Total Points: 10**<br><br>

## Advanced Exercises:

1. Use PyDESeq2 to perform a test of differential expression between participants <50 vs. >=50 years of age, controlling for sex as a covariate.
2. Use PyDESeq2 to perform a test of differential expression between participants who were (`DTHHRDY == 0`) versus were not (`DTHHRDY != 0`) on a ventillator immediately prior to death, controlling for sex and age category as covariates.<br><br>

## Additional Resources

Here are some awesome resources for you. We don't expect you to read these all, but they're relevant for discussions we had in today's lecture and could be helpful for you in your future research.

* [RNA-sequencing overview](https://www.nature.com/articles/s41576-019-0150-2). Figure 2 is especially useful and presents some of the quantification steps/tools that we haven't shown you.
* [Specific example pipeline](https://www.nature.com/articles/nprot.2016.095) from Steven Salzberg & Co, commonly used in previous iterations of bootcamp
* [Batch effects](https://www.biorxiv.org/content/10.1101/025528v1.full.pdf) discussion from Stephanie Hicks, specifically concerning single cell RNAseq
* [Replicates vs Depth](https://academic.oup.com/bioinformatics/article/30/3/301/228651) and the relation to statistical power
* The [Omnigenic inheritance](https://pubmed.ncbi.nlm.nih.gov/31051098/) model related to the discussion of cis- and trans- effects<br><br>
