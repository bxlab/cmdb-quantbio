# Bulk RNA-seq and differential expression analysis
Assignment Date: Friday, Nov. 17, 2023 <br>
Due Date: Friday, Dec. 1, 2023 <br>

## Lecture

**Slides** are available here: [Lecture slides](https://www.dropbox.com/scl/fi/qrdlsfg054xz3mkicatec/20231117_qblab_gex.pptx?rlkey=9o684hms6niwdgthanalbk81k&dl=0)


## Assignment Overview

Today, you will be examining gene expression data from the [Genotype-Tissue Expression (GTEx) Project](https://doi.org/10.1126/science.aaz1776)--a prominent resource in the field of human genetics comprising whole-genome and RNA-seq data from >800 post-mortem individuals across >50 tissues. Not all individuals donated all tissues, but there are more than 17,000 total RNA-seq samples, offering detailed insight into gene expression patterns across tissues and individuals, as well as the DNA sequence variants that modulate these differences.

One extremely common use of RNA-seq data is to compare patterns of gene expression across different conditions (i.e., predictor variables), which may include features such as sex, tissue type, genotype (e.g., wild-type versus mutant), drug treatment, time, etc. Such statistical analyses are termed "differential expression" tests, as we effectively loop over each gene and ask whether expression of that gene is correlated with the relevant predictor variable. In some study designs, you may want to "model out" or account for the effects of certain predictor variables (termed "covariates") while focusing attention on other predictors. From the statistical perspective, covariates are typically treated in the same way as other predictor variables, we simply place less focus on their effects in downstream analysis and interpretation as they may be less relevant to our key hypotheses.

Unfortunately, the most popular software packages for performing differential expression analysis (edgeR and DESeq2) were written as libraries for the R programming language, so until this year, it was not possible to introduce these packages in QuantBio Lab. As of this year, DESeq2 has been ported to Python and released as a package called "PyDESeq2". While this provides key functions of DESeq2, the package is still poorly documented, so I would like to start by motivating this exercise with a simpler (but flawed) approach for differential expression analysis based on simple linear regression. Remember that DESeq2/PyDESeq2 are also linear regression(!), but with some additional bells and whistles that I reviewed in the lecture.

## Data

The data you'll be using today is derived from the paper mentioned above, focusing on the RNA-seq samples from whole blood (755 total individuals). They were downloaded directly from the [GTEx portal](https://gtexportal.org/home/downloads/adult-gtex#bulk_tissue_expression) and slightly reformatted to save you some time on tedious data wrangling. Download from the Dropbox links below:

[Subject-level metadata](https://www.dropbox.com/scl/fi/zidlbn4rlvyv43k022mmn/gtex_metadata.txt?rlkey=j6aidakljr0739tnnzvpbg0gn&dl=0) <br>
[Gene expression matrix](https://www.dropbox.com/scl/fi/7iengpyrevd356dfq53pg/gtex_whole_blood_counts_formatted.txt?rlkey=l5h12cyher33kkzlrwi4qwf8g&dl=0)

Create a Python script for this assignment. Everything you'll need to do for this assignment will be done in this script and submitted via GitHub.


## Exercises

There are three exercises in this assignment:

1. Manually perform regression testing for differential expression
2. Perform the same tests using the `pydeseq2` library
3. Data visualization 

### Exercise 1: Perform a "homemade" test for differential expression between the sexes

#### **Step 1.1**: Loading data and importing libraries

Because the data are rectangular data of mixed type, we use Pandas to read the data into a data frame. Numpy is another alternative (and is superior in many ways), but the PyDESeq2 documentation uses Pandas, so we will stick with that for now. We will also use statsmodels to perform our "homemade" differential expression test, followed by PyDESeq2 to apply the more sophisticated test. Here is code to load those up, along with the data. You will need to modify the file paths to direct it to wherever you stored the data.

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

#### **Step 1.2**: Normalization

Before performing this test, we will first use PyDESeq2 to perform normalization across samples to account for differences in sequencing depth and RNA composition (i.e., some highly expressed genes "eating up" lots of reads and distorting the patterns). This will produce a "normalized counts" matrix that we can use for our own test. <span style="color:red;font-weight:bold">IMPORTANT NOTE:</span>: DESeq2 doesn’t actually use these normalized counts as input, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). So when you proceed to Step 3, use `counts_df`, not `counts_df_normed`!

```
counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]
```

Now extract the expression data from a given gene (let's start with the first column, DDX11L1) and merge it with the metadata, which contains the predictor variable of interest (sex). Note that sex here is encoded as 1 and 2, where 1 refers to males and 2 refers to females.

```
counts_gene = counts_df_normed.iloc[:, 0:1]
counts_gene.columns = ["counts"]
counts_gene = counts_gene.merge(metadata, on = "SUBJECT_ID")
```

#### **Step 1.3**: Statistical testing 

Use statsmodels to perform the statistical test, using `log(counts + 1)` as the response variable and sex as the predictor variable. Extract and examine the slope and p-value.

```
mod = smf.ols(formula = 'np.log(counts + 1) ~ SEX', data = counts_gene)
res = mod.fit()
slope = res.params[1]
pval = res.pvalues[1]
```

#### Step 1.4: Extend this test to all genes

Write a for-loop in Python to extend this test to all genes in your matrix. For each gene that you test, store the slopes and p-values, along with the associated gene names in one or more data structures (pandas DataFrame, arrays, ..., ?). Use the [Benjamini-Hochberg procedure](https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.fdrcorrection.html) to determine which genes are sex-differentially expressed at an FDR of 10%.

### Exercise 2: Repeat differential expression testing with PyDESeq2

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

Note that the `padj` column of `results` reports the FDR-adjusted p-value (i.e., "q-value"). The rows with a `padj` < 0.1 are the genes that are differentially expressed at an FDR of 10%. Compare these genes to those you identified in Step 1. What is the percentage of overlap? Compute this percentage in your code as a "Jaccard index", which is defined as the intersection divided by the union: `((number of genes that were significant in steps 1 and 2) / (number of genes that were significant in steps 1 or 2)) * 100%`

### Exercise 3: Visualization

Use matplotlib to create a "Volcano" plot depicting your differential expression results from Step 2. This is just another name for a scatter plot, where the x-axis shows the log2FoldChange and the y-axis shows the -log10(padj). Highlight the genes that are significant at a 10% FDR and for which the absolute value of the log2FoldChange is greater than 1 in a separate color. 

## Submission

* **Python scripts**
  * Implementation of manual DE test (Exercise 1.2-1.3) **(1 point)**
  * Code to perform DE test on all genes (Exercise 1.4) **(1 point)**
  * Running PyDESeq2 on all genes (Exercise 2) **(1 point)**
  * Code for FDR correction (Exercise 2) **(1 point)**
  * Code for percent overlap in results between methods (Exercise 2) **(1 point)**
  * Code for volcano plot (Exercise 2) **(1 point)**
* **Text submissions**
  *  List of differentially expressed transcripts (10% FDR) from Step 1 **(1 point)**
  *  List of differentially expressed transcripts (10% FDR) from Step 2 **(1 point)**
  *  Percent overlap between manual testing and PyDESeq2 **(1 point)**
* **Images**
  * Nicely formatted and labelled volcano plot **(1 point)**

## Advanced Exercises:

1. Use PyDESeq2 to perform a test of differential expression between participants <50 vs. >=50 years of age, controlling for sex as a covariate.
2. Use PyDESeq2 to perform a test of differential expression between participants who were (`DTHHRDY == 0`) versus were not (`DTHHRDY != 0`) on a ventillator immediately prior to death, controlling for sex and age category as covariates.

## Additional Resources

Here are some awesome resources for you. We don't expect you to read these all, but they're relevant for discussions we had in today's lecture and could be helpful for you in your future research.

* [RNA-sequencing overview](https://www.nature.com/articles/s41576-019-0150-2). Figure 2 is especially useful and presents some of the quantification steps/tools that we haven't shown you.
* [Specific example pipeline](https://www.nature.com/articles/nprot.2016.095) from Steven Salzberg & Co, commonly used in previous iterations of bootcamp
* [Batch effects](https://www.biorxiv.org/content/10.1101/025528v1.full.pdf) discussion from Stephanie Hicks, specifically concerning single cell RNAseq
* [Replicates vs Depth](https://academic.oup.com/bioinformatics/article/30/3/301/228651) and the relation to statistical power
* The [Omnigenic inheritance](https://pubmed.ncbi.nlm.nih.gov/31051098/) model related to the discussion of cis- and trans- effects
