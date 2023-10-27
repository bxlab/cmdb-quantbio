# Genome-Wide Association Studies

## Background

One of the central goals of human genetics is clarifying the relationship between genotypes and phenotypes. Genome-wide association studies (GWAS) emerged around 20 years ago as a useful approach for discovering genetic variation that underlies variation in human traits. At their core, GWAS involve fitting linear models to test for relationships between polymorphisms (often SNPs) and phenotypes using data from large samples of individuals.

Genotypes of nearby SNPs are correlated--a phenomenon termed "linkage disequilibrium" or "LD". This is both a blessing and a curse for GWAS. On one hand, LD means that we need not genotype every single SNP to discover associations. We merely need to genotype "tag SNPs" that segregate in LD with variants that causally influence the phenotype. On the other hand, this also means that even when we find a signficant association, it is often challenging to disentangle the causal gene and/or variant that drives the association.<br><br>

## Assignment Overview

Today we will conduct a GWAS using a materials adapted from workshop by Heather Wheeler from Loyola University Chicago (https://github.com/hwheeler01/GWAS_workshop). Note that these phenotype data are simulated. The combination of genotype and phenotype data poses a privacy risk, so real genotype and phenotype data from humans are often stored in controlled-access databases such as dbGaP.

The premise for the exercise is that you are part of a company developing two new cancer drugs called GS451 and CB1908. Some individuals in early phase trials have experienced a side effect of these drugs called lymphocytopenia (low lymphocyte counts). You are now tasked with performing a GWAS on data from lymphblastoid cell lines to search for risk factors for lymphocytopenia. The phenotype that was measured is the IC50, defined as the concentration of the drug at which 50% viability occurs.<br><br>

## Data

All the data you need for this assignment is in a zipped folder here: `~/cmdb-quantbio/assignments/lab/GWAS/extra_data/gwas_data.tar.gz`. Make a working copy of this file within your weekly homework directory.

After copying the zipped folder, you'll need to extract it with `tar -zxvf <filename.tar.gz>`. You should get 3 files:
1. genotypes.vcf
2. GS451_IC50.txt
3. CB1908_IC50.txt

The phenotype data (i.e., the IC50 for each drug) are stored in the `CB1908_IC50.txt` and `GS451_IC50.txt` files. Genotype data are stored in a VCF file (`genotypes.vcf`)—a format that you should recognize by now.

As always, if your `.gitignore` is not already set up to ignore these files (the VCF and the two phenotype files), you should update it so that they are ignored. **You should NOT be uploading any of these files.**<br><br>

## Exercises

Before you begin the exercises, create a `plotting.py` script in your submission directory. You'll be using this script to visualize 1) interesting summaries of your data and 2) your GWAS results.

Also, create a `README.md` in your submission directory. You'll be recording the commands you used for this assignment in this file, and submitting it with your plots.<br><br>

### Exercise 1: Performing PCA

#### **Step 1.1**: Compute genotype PCs

First, you want to perform PCA to compute the top genotype PCs for the samples in your data set. This is useful not only for visualizing the population structure in your data set, but also allows you to control for this structure when you perform your GWAS later.

Using `plink`, perform PCA on the genotypes of the cell lines in the data set. Output the top 10 principal components.

Record the `plink` command(s) you used in your `README.md`.<br><br>

#### **Step 1.2**: Plot genotype PCs

Now that you've computed the top genotyping PCs, you can visualize the population structure in your data set.

In your `plotting.py` script, plot the top **2** genotype principal components from `plink` to visualize the genetic relatedness between the cell lines in the data set. Ensure that your plot is nicely formatted and properly labeled.<br><br>

### Exercise 2: The allele frequency spectrum

#### **Step 2.1**: Compute allele frequencies

The shape of the allele frequency specturm (AFS) can tell you interesting things about your sample, such as demographic history, and evidence of selection (read more [here](https://en.wikipedia.org/wiki/Allele_frequency_spectrum). While you will not be running these specific tests today, it is still interesting (and good practice) to visualize the AFS.

Usually, the AF of each variant is stored in the `INFO` field in your VCF. The vcf file we're using does not have an `INFO` field, so you'll need to calculate AF yourself.

Using `plink`, calculate allele frequencies for the variants in your VCF and output these to a file.

Record the `plink` command(s) you used in your `README.md`.<br><br>

#### **Step 2.2**: Plot AFS

Now that you've computed allele frequencies, you can visualize the AFS of the samples in your data set.

In your `plotting.py` script, plot the AFS of the samples in your data set as a *histogram*. Ensure that your plot is nicely formatted and properly labeled.<br><br>

### Exercise 3: GWAS

#### **Step 3.1**: Running the GWAS

Now you're ready to run your GWAS!

Using `plink`, perform quantitative association testing for each of your two phenotypes (you will be running `plink` one each per phenotype). Use the top 10 principal components you computed in **Step 1.1** as covariates in your analysis to control for non-independence due to relatedness/population structure.
      <li> Be sure to use the <code>--allow-no-sex</code> option</li>
      <li> You may find this portion of the <code>plink</code> <a href="https://zzz.bwh.harvard.edu/plink/anal.shtml">documentation</a> helpful  for performing association testing on each of the phenotypes.</li>
      <li> HINT: We do NOT have case/control phenotpe data, so `--assoc` is not the correct `plink` flag.

We encourage you to try this on your own first, but you can use the hint below to check your answer (this WILL give you the full `plink` command, be warned).

<details><summary><b>CLICK HERE IF YOU WANT TO CHECK YOUR ANSWER:</b></summary>
  <code>plink --vcf genotypes.vcf --linear --pheno &lt;<phenotype>.txt&gt; --covar &lt;pca.eigenvec&gt; --allow-no-sex --out &lt;<phenotype>_gwas_results&gt;</code>
</details>

You should have two output files, one for each phenotype.

Record the `plink` commands you used in your `README.md`.<br><br>

#### **Step 3.2**: Visualing GWAS results

Now that you've run your GWASes, you want to visualize the results.

In your `plotting.py` script, produce a two-panel figure (one column, two rows) depicting the Manhattan plot for each of your two GWAS analyses (i.e. each phenotype). Each Manhattan plot should be one of the two panels in your figure. In each Manhattan plot, highlight SNPs with p-values less than 10<sup>-5</sup> in a different color. 

Ensure that your figure/panels are nicely formatted and properly labeled.<br><br>

#### **Step 3.3**: Visualizing effect-size

You want to dig deeper into the top GWAS hits from your analyses.

Choose one of the traits for which you performed GWAS.

In your `plotting.py` script: for the top associated SNP (lowest p-value) of that trait, plot the effect size of that variant on the chosen trait by creating a boxplot of the phenotype stratified by genotype.<br><br>

#### **Step 3.4**: What gene could it be?

For the top loci associated with **each** of your two phenotypes, use the <a href="http://genome.ucsc.edu/cgi-bin/hgGateway">UCSC Genome Browser</a> to investigate the potential causal genes in the region. Note that the reference genome build being used here is **hg18**. (How could you figure this out if you didn't know?).

In your `README.md`, summarize what you found (i.e. what gene(s) are closest to each top hit? How might these genes impact the trait?).<br><br>

<!-- ### Hints

- You can perform the PCA and allele frequency calculations in Python, but `plink` can perform both as well. Either way, you will probably want to use matplotlib to produce the plots.
- If you do PCA in python, be sure you have one datapoint per individual, as opposed to per SNP.
- Each SNP only has one allele frequency.
- Notice how the IDs are encoded in the phenotype files [phenotypes](https://www.cog-genomics.org/plink2/input#pheno]). `plink` expects both a family ID and a sample ID. In your VCF, they are separated by an underscore, but in your phenotype file they are separated by a tab-character. The files we provided are already formatted in this way, but keep this in mind for future work where reformatting may be necessary. -->

## Submission

For this assignment you should submit the following:

1. `plink` code to perform PCA and top 10 PCs (**1pt**)
2. PCA figure (**1pt**)
3. Code to generate allele frequency histogram (**1pt**)
4. Allele frequency histogram + clear labels (**1pt**)
5. Script to perform association testing using `plink` (**1pt**)
6. Manhattan Plots (**1pt each**)
7. Phenotype effect boxplot (**1pt**)
8. Written summary about potential causal variants in `README.md` (**1pt**)

**DO NOT** push any raw data to Github unless we explicitly ask for it!

## Advanced Exercises

1. In addition to your Manhattan plot, produce a QQ plot for each of the two traits, summarizing the results.

2. Examine the frequency of the top associated variants. What is its frequency in various global populations? Use the GGV Browser to explore this: https://popgen.uchicago.edu/ggv

3. Is the variant itself necessarily causal in driving the association? Investigate the local haplotype structure using the LDLink website's "LDproxy" tool: https://ldlink.nci.nih.gov/?tab=ldproxy




