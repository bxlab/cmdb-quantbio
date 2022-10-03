# Assignment 4: Genome-Wide Association Studies
Assignment Date: Friday, Sept. 30, 2022 <br>
Due Date: Friday, Oct. 7, 2022 @ 1:00pm ET <br>

**Please do a `git pull` within your `~/cmdb-quantbio/` directory**

## Lecture

Slides are available here: [Lecture Slides](https://github.com/bxlab/cmdb-quantbio/blob/main/assignments/lab/GWAS/slides_asynchronous_or_livecoding_resources/20220930_qblab_gwas.pptx?raw=true)

## Background

One of the central goals of human genetics is clarifying the relationship between genotypes and phenotypes. Genome-wide association studies (GWAS) emerged around 20 years ago as a useful approach for discovering genetic variation that underlies variation in human traits. At their core, GWAS involve fitting linear models to test for relationships between polymorphisms (often SNPs) and phenotypes using data from large samples of individuals.

Genotypes of nearby SNPs are correlated--a phenomenon termed "linkage disequilibrium" or "LD". This is both a blessing and a curse for GWAS. On one hand, LD means that we need not genotype every SNP to discover associations. We merely need to genotype "tag SNPs" that segregate in LD with variants that causally influence the phenotype. On the other hand, this also means that even when we find a signficant association, it is often challenging to disentangle the causal gene and/or variant that drives the association..

## Assignment Overview

Today we will conduct a GWAS using a materials adapted from workshop by Heather Wheeler from Loyola University Chicago (https://github.com/hwheeler01/GWAS_workshop). Note that these phenotype data are simulated. The combination of genotype and phenotype data poses a privacy risk, so real genotype and phenotype data from humans are often stored in controlled-access databases such as dbGaP.

The premise for the exercise is that you are part of a company developing two new cancer drugs called GS451 and CB1908. Some individuals in early phase trials have experienced a side effect of these drugs called lymphocytopenia (low lymphocyte counts). You are now tasked with performing a GWAS on data from lymphblastoid cell lines to search for risk factors for lymphocytopenia. The phenotype that was measured is the IC50, defined as the concentration of the drug at which 50% viability occurs.

## Data

All the data you need for this assignment is in a zipped folder here: `~/cmdb-quantbio/assignments/lab/GWAS/extra_data/gwas_data.tar.gz`. Copy this file to the `answers` directory you made for this assignment.

After copying the zipped folder, you'll need to extract it with `tar -zxvf <filename.tar.gz>`. You should get 3 files:
1. genotypes.vcf
2. GS451_IC50.txt
3. CB1908_IC50.txt

The phenotype data (i.e., the IC50 for each drug) are stored in the `CB1908_IC50.txt` and `GS451_IC50.txt` files. Genotype data are stored in a VCF file (`genotypes.vcf`)—a format that you should recognize by now.

You shouldn't be pushing any of these files to your remote repo, so update your `.gitignore` file as necessary.

## Assignment

1. Create a `README.md` file where you'll record all of the commands you used for this assignment.
2. Using `plink`, perform PCA on the genotypes of the cell lines. **Output at least the top 10 principal components.** Then, visualize genetic relatedness between the cell lines by plotting the first two components. You'll need to submit this figure, as well as the commands you used to generate the principal components in your `README.md`.
3. Visualize the allele frequency spectrum by plotting a histogram of allele frequencies. You'll need to submit this figure, as well as the commands you used to generate the AFs in your `README.md`.
  - The vcf file we're using does not have an `INFO` field, so you'll need to calculate AF yourself. You can do this in Python, but `plink` can perform this calculation as well, which may be easier. Your choice.
4. Using `plink`, perform quantitative association testing for each phenotype. Use the top 10 principal components (eigenvectors) as covariates in your analysis, to adjust for non-independence due to relatedness. Record the commands you used in your `README.md`
  - Be sure to use the `--allow-no-sex` option
  - You may find this portion of the [`plink` documentation](https://zzz.bwh.harvard.edu/plink/anal.shtml) helpful for performing association testing on each of the phenotypes.
***
<details><summary> HINT: </summary>
`plink --vcf genotypes.vcf --linear --pheno <phenotype.txt> --covar <pca.eigenvec> --allow-no-sex --out <phenotype_gwas_results>`
</details>
***
5. For each phenotype, produce a Manhattan plot. In the Manhattan plot, highlight SNPs with p-values less than 10<sup>-5</sup> in a different color. You'll need to submit this figure.
6. Choose one of the traits for which you performed GWAS. For the top associated SNP, visualize the effect size by creating a boxplot of the phenotype stratified by genotype.
7. For the top loci associated with each of the phenotypes, use the UCSC Genome Browser (http://genome.ucsc.edu/cgi-bin/hgGateway) to investigate the potential causal genes in the region. Note that the reference genome build being used here is hg18. (How could you figure this out if you didn't know?). Summarize your results in your `README.md`.

### Hints

- You can perform the PCA and allele frequency calculations in Python, but `plink` can perform both as well. Either way, you will probably want to use matplotlib to produce the plots.
- If you do PCA in python, be sure you have one datapoint per individual, as opposed to per SNP.
- Each SNP only has one allele frequency.
- Notice how the IDs are encoded in the phenotype files [phenotypes](https://www.cog-genomics.org/plink2/input#pheno]). `plink` expects both a family ID and a sample ID. In your VCF, they are separated by an underscore, but in your phenotype file they are separated by a tab-character. The files we provided are already formatted in this way, but keep this in mind for future work where reformatting may be necessary.


## Submit

For this assignment you should submit the following:

1. Your `README.md` file containing the commands you used, as well as your answer for Part 7
2. Any python scripts you wrote
3. The plots for Parts 2, 3, 5, and 6

**DO NOT** push any raw data to Github unless we explicitly ask for it!

## Advanced Exercises

1. In addition to your Manhattan plot, produce a QQ plot for each of the two traits, summarizing the results.

2. Examine the frequency of the top associated variants. What is its frequency in various global populations? Use the GGV Browser to explore this: https://popgen.uchicago.edu/ggv

3. Is the variant itself necessarily causal in driving the association? Investigate the local haplotype structure using the LDLink website's "LDproxy" tool: https://ldlink.nci.nih.gov/?tab=ldproxy




