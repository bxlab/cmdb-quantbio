## Examining the causes of allele specific expression

Genetic variation near a gene (e.g., within the promoter) may contribute to differences in gene expression among individuals. GTEx has used different approaches to assess evidence of such *cis*-regulatory effects on gene expression. Two complementary approaches for discovering these effects are expression quantitative trait locus (eQTL) mapping and allele-specific expression (ASE) analysis.

In eQTL mapping, subject genotypes at single nucleotide polymorphisms (SNP) near a gene are correlated with the expression level of that gene (mRNA abundance, measured by the normalized count of aligned reads) across a population sample. Genes with at least one detected correlated SNP (eQTL) are called eGenes.

In contrast, ASE typically focuses on individual samples and asks whether the maternally inherited haplotype (which may be tagged by one or more heterozygous SNPs) is expressed at a different level than the paternally inherited haplotype. While the identity of the maternal versus paternal haplotype is typically not known unless the parents are also sequenced, a difference in expression level between the observed alleles of one or more SNPs that compose a haplotype suggests the existence of a linked regulatory variant in *cis*.

GTEx recently published some of the first [long-read sequencing data](https://gtexportal.org/home/downloads/adult-gtex/long_read_data) from human samples and quantified expression of maternally versus paternally inherited haplotypes with a new tool called LORALS.

One question that we might ask is whether genes that exhibit evidence of ASE based on the long-read sequencing data are enriched for eQTLs that were previously detected with short-read sequencing data. We will address this question in the following exercise.

Please submit your answers as a single `.R` script with comments that separate out answers to each question below. (Remember that comment lines start with `#` and are ignored by R). For questions regarding interpretation or discussion of your results, please include your answers as comments interspersed with your code.

1.  Load the `tidyverse` and `qvalue` packages, and use the function `read_delim()` to read in the long-read-based allele-specific expression data (`GTEx_lorals_ase_data.txt`).


2.  View the first rows of the tibble by simply entering the variable name in which you stored it. Explore the data a bit. How many unique samples are represented? How many unique genes? How many unique tissues?

    Note that the column `is_eGene` denotes whether the gene possessed a *cis*-eQTL based on short-read RNA-seq data from GTEx from the same tissue (cultured fibroblasts). What are the counts of `TRUE` and `FALSE` in this column? Ignore the fact that some genes are represented multiple times, as their allele specific expression is measured in multiple samples.


3.  The column `h1` contains the counts of long-reads that were inferred to originate from one haplotype (e.g., the maternally-inherited haplotype), whereas the column `h2` contains the counts of long-reads that were inferred to originate from the other haplotype. Reads were assigned to haplotypes based on observation of alleles of one or more heterozygous exonic SNPs.
    -   What is the minimum count of `h1`? What is the maximum? What is the minimum count of `h2`? What is the maximum?

    -   Create a scatter plot of `h1` versus `h2`. Use `geom_abline()` to add a diagonal line indicating x = y (i.e., perfect allele balance).

    -   Re-plot the data while scaling the x and y axes using `scale_x_log10()` and `scale_y_log10()`.


4.  The following function applies a binomial test and returns a p-value, quantifying evidence of allele specific expression for a given gene in a given sample (i.e., a row of your data). Apply this function to each row of your tibble using the `rowwise() %>% mutate()` chain of functions, adding a column that contains the p-values. Read about how to use these functions together [here](https://dplyr.tidyverse.org/articles/rowwise.html).\

    ``` R
    ase_pval <- function(hap1_counts, hap2_counts) {   
      return(binom.test(hap1_counts, (hap1_counts + hap2_counts), p = 0.5)$p.value) 
    }
    ```

5.  Copy you log-scaled scatter plot from question 3, but now color the points based on their p-value.


6.  Use the `qvalue()` function to convert your p-values to q-values. Note that q-values are stored in the `$qvalue` slot of the object output by this function.

    The q-value is defined as the minimum false discovery rate (FDR) at which an observation is deemed significant. How many of your gene-sample pairs (i.e., rows) exhibit significant evidence of ASE at a 10% FDR? This means that we expect our list of significant hits to contain 10% false positives (though we don't know which particular observations are false positives.


7.  We are now ready to test whether our genes with evidence of ASE are enriched for eGenes (i.e., genes with eQTLs). We will use permutation to perform this test.

    How many of your observations with significant ASE (at a 10% FDR) based on your analysis above are also eGenes? Hint: the function `nrow()` can count the number of rows in a data frame or tibble.


8.  To figure out if the value you computed above represents a significant association between ASE and eGenes, we need to construct a null distribution, which we will accomplish by permutation (i.e., shuffling eGenes with respect to the other columns of your data frame).

    Use the `mutate()` function with the `sample()` fuction to randomize the order of `is_eGene` with respect to the rest of your data frame and re-compute your answer to question 7. Use a for loop to repeat this procedure 1000 times to produce a null distribution. Plot this null distribution as a histogram.


9.  Use `geom_vline()` to overlay your observed value from question 7 on the histogram from question 8. \

    Compute a one-sided p-value by asking what proportion of your null simulations are as extreme or more extreme than your observed value. What is your statistical and biological interpretation of this result?

