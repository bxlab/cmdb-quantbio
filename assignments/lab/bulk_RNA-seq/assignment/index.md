# RNA-seq differential expression analysis

## Assignment Overview

Today, you will be examining gene expression data from the [Genotype-Tissue Expression (GTEx) Project](https://doi.org/10.1126/science.aaz1776)--a prominent resource in the field of human genetics comprising whole-genome and RNA-seq data from >800 post-mortem individuals across >50 tissues. Not all individuals donated all tissues, but there are more than 17,000 total RNA-seq samples, offering detailed insight into gene expression patterns across tissues and individuals, as well as the DNA sequence variants that modulate these differences.

One extremely common use of RNA-seq data is to compare patterns of gene expression across different conditions (i.e., predictor variables), which may include features such as sex, tissue type, genotype (e.g., wild-type versus mutant), drug treatment, time, etc. Such statistical analyses are termed "differential expression" tests, as we effectively loop over each gene and ask whether expression of that gene is correlated with the relevant predictor variable. In some study designs, you may want to "model out" or account for the effects of certain predictor variables (termed "covariates") while focusing attention on other predictors. From the statistical perspective, covariates are typically treated in the same way as other predictor variables. We simply place less focus on their effects in downstream analysis and interpretation as they may be less relevant to our key hypotheses.

The most popular software packages for performing differential expression analysis (edgeR and DESeq2) are written in R, so we will be working in R for this assignment. Before running this package, however, I would like to start by motivating this exercise with a simpler (but somewhat flawed) approach for differential expression analysis based on simple linear regression. Remember that DESeq2 is also linear regression(!), but with some additional bells and whistles that I reviewed in the lecture.<br><br>

## Data

The data you'll be using today is derived from the GTEx Project mentioned above, focusing on the RNA-seq samples from whole blood. They were downloaded directly from the [GTEx portal](https://gtexportal.org/home/downloads/adult-gtex#bulk_tissue_expression) and slightly reformatted to save you some time on tedious data wrangling. I also downsampled the data to only 106 individuals to ensure that the analysis runs fast on your laptops.


Download from the Dropbox links below:

[Subject-level metadata](https://www.dropbox.com/scl/fi/vz0b6ybc66uqlg8c1pv0u/gtex_metadata_downsample.txt?rlkey=gw88lmucf7xc0z8yh3p0qb2rl&dl=0)

[Gene expression data](https://www.dropbox.com/scl/fi/aryv6uofzn329o42e7dst/gtex_whole_blood_counts_downsample.txt?rlkey=ojh88vqelrck9j2v0h2alfv6k&dl=0)

[Locations of genes on chromosomes](https://www.dropbox.com/scl/fi/2qz19ctbubu3d0xadn20e/gene_locations.txt?rlkey=546b2clj8jnvgvvr7wzazxgvk&dl=0)<br><br>

## Exercises

There are three exercises in this assignment:

1. Load the data and visualize top principal components to explore sample similarity
2. Perform differential expression analysis using the `DESeq2` library
3. Visualize and interpret your results

Before you do anything else, create an R script for this assignment. Everything you'll need to do for this assignment will be done in this script and submitted via GitHub. Use comments to separate your code blocks into Exercise 1, Step 1.1; Exercise 1, Step 1.2; etc., and include written answers to questions as comments.<br><br>

### Exercise 1: Perform principal component analysis

A reasonable initial step for many gene expression analyses is to apply an unsupervised method to examine which samples are more or less similar to one another in their broad patterns of gene expression. We will implement this step using functions built into the `DESeq2` package.

#### **Step 1.1**: Loading data and importing libraries

**1.1.1**: First, load the `DESeq2` and `tidyverse` libraries, and use the `read_delim()` function to read the data and metadata into separate tibbles (the `tidyverse` version of data frames). 

**1.1.2**: The gene expression file reports individuals (subjects) in columns and genes in rows. Each entry in the table records the number of reads mapping to that particular gene for the whole-blood sample from that particular subject. Use the function `column_to_rownames()` to move the gene names, which are stored in the first column to instead be stored as rownames. This is a necessary step due to the format that DESeq2 expects the data, whereby the formal entries in the tibble should be solely numeric.

**1.1.3**: The metadata table records the sex, age, and cause of death for each sample. Ages are encoded as integers, ranging from 0 to 4 that represent decade-long age intervals, where 0 denotes ages 20-29, 1 denotes ages 30-39, ..., and 5 denotes ages 70-79. By keeping encoding this variable as an integer in our model, rather than a categorical variable, we are modeling linear, monotonic relationships between gene expression level and age. Use the function `column_to_rownames()` to move the subject IDs, which are stored in the first column to instead be stored as rownames. 

**1.1.4**: After loading and reformatting as instructed above, peek at the first few rows and columns of the tibbles, just to understand how they are organized and ensure that they were loaded correctly.

<br>

#### **Step 1.2**: Create a DESeq2 object

**1.2.1**: It is critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. Fortunately, this is already the case for the data we have provided, but it is important to double-check! Write a line or two of code to perform this check.

**1.2.2**: Next, load the data into a `DESeq2` object using the `DESeqDataSetFromMatrix()` function. Name this object `dds`. This function takes three required arguments: 
- `countData`: the gene expression count matrix (or tibble)
- `colData`: the metadata tibble
- `design`: the model formula, which uses the linear model notation that we learned from bootcamp, but you can leave off the response variable, since that will always just be the expression of a given gene. For example, if my colData had a variable of interest called CONDITION and another covariate called SMOKING_STATUS, I might want to specify my formula as `design = ~ CONDITION + SMOKING_STATUS`.

Replace the brackets below with your code. For your model formula, include the variables SEX, DTHHRDY, and AGE. (The brackets themselves can be removed).

```r
# load DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = <>,
                              colData = <>,
                              design = <>)
```

#### **Step 1.3**: Normalization and PCA

**1.3.1**: As we did last week, we will want to apply a variance stabilizing transformation (VST) before performing principal component analysis, which has the effect of giving genes with different total expression levels more equal weights. Use the `vst()` function to apply VST to your DESeq2 object, storing it in a new variable called `vsd` (for variance-stabilized data).

**1.3.2**: Now apply the `plotPCA()` function to your VST-normalized data. Perform this step three separate times, each time providing a different value to the `intrgroup =` argument to color the points by various variables in your metadata. Save your plots as PNG files and include them in your submission.

**1.3.3**: What proportion of variance in the gene expression data is explained by each of the first two principal components? Which principal components appear to be associated with which subject-level variables? Interpret these patterns in your own words and record your answers as a comment in your code.

### Exercise 2: Perform differential expression analysis

### Step 2.1: Perform a "homemade" test for differential expression between the sexes

You're curious whether any genes are differentially expressed between sexes (in whole blood). Before you do the DE analysis the "right" way with DESeq2, you're first going to perform your own naive "homemade" analysis to get an idea of what DESeq2 is doing behind the scenes at a basic level. This will involve regressing expression quantifications from GTEx onto sample sex, and identifying genes with significant differences in expression between sexes. We will use the VST-normalized data for this test, but note that in our later analyses, we will go back to the raw counts!

As a first step, we can use the code below to extract the VST expression matrix and bind it to the metadata, such that the whole data frame can be used as input to a regression model:

```r
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()

vsd_df <- bind_cols(metadata_df, vsd_df)
```

For example, we can test for differential expression of the gene WASH7P by running the following code:

```r
m1 <- lm(formula = WASH7P ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()
```

**2.1.1**: Does WASH7P show significant evidence of sex-differential expression (and if so, in which direction)? Explain your answer.

**2.1.2**: Now repeat your analysis for the gene SLC25A47. Does this gene show evidence of sex-differential expression (and if so, in which direction)? Explain your answer.

#### **Step 2.2**: Perform differential expression analysis "the right way" with DESeq2

**2.2.1**: Now apply the `DESeq()` function to fit the regression model that you previously specified in your `design =` argument to each gene in the dataset. This function is a wrapper function that actually performs several different steps. Note that the function should be run on the original un-normalized count data (`dds`) rather than the data that underwent VST. This is because `DESeq2()` is designed to model the counts themselves and already includes a step to estimate library sizes (i.e., the total amount of data per sample) and correct for this in the model. Store the results back in the same object, `dds`.

#### **Step 2.3**: Extract and interpret the results for sex differential expression

Results of this differential expression analysis can be extracted using the `results()` function. By default, this function will focus on whatever variable was specified last in your `design =` formula, but while controlling for all of the other covariates listed in the model formula. You can override the default behavior and tell `results()` to focus on any variable of interest that was included in that formula. For example, if my original model formula was `design = ~ CONDITION + SMOKING_STATUS`, and `CONDITION` had the levels `treatment` and `control`, I could extract differential expression results for `CONDITION` (controlling for `SMOKING_STATUS`) and store it in a tibble using the following code.

```r
res <- results(dds, name = "CONDITION_treatment_vs_control")  %>%
  as_tibble(rownames = "GENE_NAME")
```

**2.3.1**: Extract the differential expression results for the variable SEX.

**2.3.2**: The `padj` column of `results` reports the FDR-adjusted p-value (i.e., "q-value"). The rows with a `padj` < 0.1 are the genes that are differentially expressed at an FDR of 10%. A lot of rows will have missing `padj` values. This is because DESeq2 does a filtering step to remove genes that it thinks have a low probability of being DE, and only performs FDR correction for the set of genes it does not filter out. You can ignore genes with missing `padj` values.

How many genes exhibit significant differential expression between males and females at a 10% FDR?

**2.3.3**: Now load the mappings of genes to chromosomes, which have already downloaded from Dropbox (see top of page). Merge these data with your differential expression results using the `left_join()` function, setting `by = GENE_ID` (i.e., the shared column on which you want to join). Order the merged tibble by `padj`, from smallest to largest.

Examine your top hits. Which chromosomes encode the genes that are most strongly upregulated in males versus females, respectively? Are there more male-upregulated genes or female-upregulated genes near the top of the list? Interpret these results in your own words. 

**2.3.4**: Examine the results for the two genes (WASH7P and MAP7D2) that you had previously tested with the basic linear regression model in step 2.1. Are the results broadly consistent?

#### **Step 2.4**: Extract and interpret the results for differential expression by death classification

**2.4.1**: Now repeat your analysis from step 2.2 above, but focusing on death classification as your variable of interest. This can be accomplished with the `results()` function, specifying the argument `name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes"`.

How many genes are differentially expressed according to death classification at a 10% FDR?

**2.4.2**: Interpret this result in your own words. Given your previous analyses, does it make sense that there would be more genes differentially expressed based on type of death compared to the number of genes differentially expressed according to sex?

### Exercise 3: Visualization

**3.1**: Use `ggplot2` to create a "Volcano" plot depicting your differential expression results from Exercise 2.3. A volcano plot is a scatter plot, where the x-axis shows the `log2FoldChange` and the y-axis shows the -log10(`padj`).

Highlight the genes that are significant at a 10% FDR **AND** for which the absolute value of the log2FoldChange is greater than 1 in a separate color.

Output this plot to a `.png` that you will upload with your assignment.<br><br>

## Submission

1. R script to run DE analysis **(8 points total)**
  * Code to load data and create a DESeq2 object (Steps 1.1 and 1.2) **(1 point)**
  * Code to perform VST and PCA (Step 1.3) **(1 point)**
  * Interpretation of PCA (Step 1.3 questions) **(1 point)**
  * Code to apply differential expression analysis (Step 2.1) **(1 point)**
  * Code to extract results for sex differential expression (Step 2.2) **(1 point)**
  * Interpretation of sex differential expression (Step 2.2) **(1 point)**
  * Code to extract results for differential expression by death classification (Step 2.3) **(1 point)**
  * Interpretation of sex differential expression by death classification (Step 2.4) **(1 point)**

2. Plots **(2 points total)**
  * Nicely formatted and labelled PCA plots from Exercise 1 **(1 points)**
  * Nicely formatted and labelled volcano plot from Exercise 2 **(1 points)**

**Total Points: 10**<br><br>

## Advanced Exercises:

1. Search the literature to interpret some of your results from the differential expression analysis. Aside from the genes on sex chromosomes, do other genes that you see turning up as strongly differentially expressed make sense in light of the literature?

2. Use DESeq2 to perform a test of differential expression with respect to age. How many genes are differentially expressed? How do you interpret genes with positive versus negative associations?

3. Try applying functional enrichment analysis to ask if your genes are overrepresented for some particular functional category. The packages `fgsea` and `msigdbr` provide one nice approach (see tutorial here: https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html). 

## Additional Resources

Here are some awesome resources for you. We don't expect you to read these all, but they're relevant for discussions we had in today's lecture and could be helpful for you in your future research.

* [RNA-sequencing overview](https://www.nature.com/articles/s41576-019-0150-2). Figure 2 is especially useful and presents some of the quantification steps/tools that we haven't shown you.
* [Specific example pipeline](https://www.nature.com/articles/nprot.2016.095) from Steven Salzberg & Co, commonly used in previous iterations of bootcamp
* [Batch effects](https://www.biorxiv.org/content/10.1101/025528v1.full.pdf) discussion from Stephanie Hicks, specifically concerning single cell RNAseq
* [Replicates vs Depth](https://academic.oup.com/bioinformatics/article/30/3/301/228651) and the relation to statistical power
* The [Omnigenic inheritance](https://pubmed.ncbi.nlm.nih.gov/31051098/) model related to the discussion of cis- and trans- effects<br><br>
