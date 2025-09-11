# Genome Content

## Assignment Overview

The goal of today's lab is to learn about what types of information exists within a genome, the formats that data can take, how to visualize and obtain publicly available data, and how to use various data sources to identify relevant subsets of a genome. In order to carry out today's assignment, you will be making extensive use of the [UCSC Genome Browser](https://genome.ucsc.edu/) for data visualization, filtering and downloading data, and analysis. The aim of this assignment is to look at the prevalence (and therefor the infered genomic tolerance for) SNPs within various genomic features and across varying SNP frequencies in the population.

**Important**
Please commit and push as you finish each part (in the case of 2.1, commit and push as you do each task within the step). This will really help us gauge the where people are being challenged and better help you all and future students.


## Data

 You will be provided with a set of bed files containing a set of common SNPs for human chromosome 1 (genome build hg38) partitioned by minor allele fequency (MAF). These SNPs were obtained from the UCSC Genome Browser from the dbSNP release 151 common SNPs. There is also a bed file dividing chromosome 1 into contiguous 20Kb bins.

## Exercises

### Exercise 1: Obtain the data

#### **Step 1.1** Download data

You'll start by downloading the SNP files from [here](https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/assignments/lab/genome_content/extra_data/chr1_snps.tar.gz). You will need to unpack the contents using the following command:

```bash
tar xzf chr1_snps.tar.gz
```

You should now see 6 files, chr1_snps_0.*.bed for values 1-5 and genome_chr1.bed.


#### **Step 1.2** Download feature files using UCSC Genome Table Browser

Using the Table Browser tool in the genome browser, download the GenCode v46 known genes from the `Genes and Gene Predictions` group. Make sure you restrict the region of interest to chromosome 1. You will also want to specify under the `Output format` menu that you only want `Selected fields from the primary and related tables`. Also, specify a file name to download the data to. Then press the `Get output` button. You can now select the specific fields you want. For this assignment, you only want the chromosome, txStart (transcription start), and txEnd (transcription end) positions.

Repeat the above process, except this time select chromosome, exonStarts, and exonEnds. This will give you a bed file with all known coding exons.


#### **Step 1.3** Use Table Browser to obtain Encode cCREs

You will also need to download a set of regulatory elements using the Table Browser. The procedure should be identical to step 1.2 except you will find these data in the `Regulation` group under the track `Encode cCREs`. Again, you only need the chromosome, chromStart, and chromEnd fields.

#### **Step 1.4** Use bedtools to merge elements within each feature file

Unfortunately, because there are many isoforms of genes, many regions in the genome are represented multiple times in the gene and exon bed files that you just downloaded. In order to create a file without overlapping ranges, you will need to use `bedtools merge` to eliminate the redundancy in each file. If you simply try to run bedtools on each file as is, you will get an error telling you that the file needs to be sorted by chromosome and coordinate. You can do this with `bedtools sort`.

For the TAs' sake, name your final feature files `<feature>_chr1.bed` where `<feature>` is each feature name (genes, exons, or cCREs) for a total of 3 files.


#### **Step 1.5** Use bedtools to create intron feature file

It would be nice to look at introns as well as exons. Because of the multiple isoform issue, many introns contain exons that are skipped in a particular isoform. To get around this, you will be using the merged gene and merged exon files that you created in step 1.4 to infer regions of the chromosome that are purely intronic by subtracting exonic intervals from the gene intervals using `bedtools subtract` to create the intron_chr1 bed file.

For the TAs' sake, name your intron feature file `introns_chr1.bed`.


#### **step 1.6** Use bedtools to find intervals not covered by other features

Finally, we want to see how SNP density changes in regions of chromosome 1 that are not covered by the other features of interest. To create an `other` interval bed file, you can use the same approach as in step 1.5, but using the `genome_chr1.bed` file as the target and subtracting the exon, intron, and cCRE files (this can be done in a single call of `bedtools subtract`). You should now have 4 feature bed files for the downstream analysis, exons, introns, cCREs, and other.

For the TAs' sake, name your intron feature file `other_chr1.bed`.


### Exercise 2: Count feature SNPs and determine enrichment

Using the datafiles you created in part 1 along with the SNP bed files partitioned by MAF, you will be finding how many SNPs overlap each, finding the mean SNP/bp for each feature and then calculating the enrichment of each feature SNP density for each MAF level.


#### **Step 2.1** Create a bash script to find the overlap of SNPs and features and calculate the SNP density enrichment for each MAF-feature combination

Write a bash script that loops through each MAF file and each feature bed file (this sounds a lot like a nested `for` loop situation) and use `bedtools coverage` to find how many SNPs fall within each set of features. You don't actually need the set of SNPs that overlap each feature set, but instead are interested in the total count of SNPs. The resulting files should have all the information you need to find the SNP density enrichments.

You will need to sum two values, the number of SNPs and the total size of the feature. Double check what the output columns of `bed coverage` are. In order to sum the values in a colum, remember that you can use the `awk` command covered in lecture, subbing in the correct column number for `$1` (i.e. column 3 would be `$3`) and putting in your file name for `<filename>`.

```bash
SUM=$(awk '{s+=$1}END{print s}' <filename>)
```

You will need to calculate the ratio of these sums to determine the SNP density. You can use `bc` for this. Just remember to use the `-l` flag to give you a decimal answer. You will also need the `-e` flag before your math equation (remember to put it in quotation marks). Also, don't forget that you want an enrichment value, not a ratio. How do you get an enrichment? Divide by the background (background = total # SNPs for given MAF / chromosome length).

Finally, you need to save the results in a file. I suggest creating the file by putting in an informative header, like:

```bash
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt
```

You should have 5 SNP files (different MAF levels) and 4 feature files, for a total of 20 combinations. Each line should have 3 elements, the MAF, the feature, and the enrichment. Please name this file `snp_counts.txt` for ease of grading.

To help you, here is a pseudo code to scaffold the general approach:

```
Create the results file with a header
Loop through each possible MAF value
  Use the MAF value to get the file name for the SNP MAF file
  Find the SNP coverage of the whole chromosome
  Sum SNPs from coverage
  Sum total bases from coverage
  Calculate the background
  Loop through each feature name
    Use the feature value to get the file name for the feature file
    Find the SNP coverage of the current feature
    Sum SNPs from coverage
    Sum total bases from coverage
    Calculate enrichment
    Save result to results file
  ```


#### **Step 2.2** Use R to plot data 

Once you have your SNP enrichment values for each MAF and feature, write an R script that loads them into a dataframe and plots them using ggplots2. The plot should have 4 lines, one for each feature with the MAF on the X-axis and the log2-transformed enrichment on the Y-axis. You should include a legend and label your axes.

Please name this plot `snp_enrichments.pdf`.


#### **Step 2.3** Draw conclusions from plot

Based on the plotted data from step 2.2, you should be able to draw some conclusions. Put answers in the file `README.md`

*1. Which feature appears to be under the highest level of purifying (negative) selection? How can you tell?*
*2. Why do you think most features start out more enriched at lower MAFs and become more depleted as the MAF increases?*
*3. Does the relative levels of enrichment/depletion of features with respect to each other make sense? Explain.*


## Advanced Excercises (Optional)

### Exercise 3: Loading custom tracks

Another thing that you can do with the genome browser is visualize your own data alongside the large number of publicly available tracks. You can also perform some analyses right in the browser.

#### **Step 3.1** Use bedtools to create a SNP coverage bedgraph file

With `bedtools coverage`, create a bedgraph file with the number of SNPs in each interval of `genome_chr1.bed` file. You will need to first concatenate all of the MAF SNP files into a single file prior to running the coverage tool. You will need to select only the relevant columns from the bedtools output using the command line tool `cut` and selecting the chromosome, start position, end position, and count. In order for the genome browser to properly interpret your bedgraph, add the following line of text to the beginning of your bedgraph file:

`track type=bedGraph name="20Kb binned SNPs" visibility="full" windowingFunction="mean"`


#### **Step 3.2** Upload your custom track

Load your new track into the genome browser. Click the `Manage custom tracks` button under the track display. Then use `Choose File` to select your bedgraph file and finally click `Submit`. You should now be able to go back to the browser and see your SNP coverage track.

Let's compare it to another track, `Cons 30 Primate` under the `Comparative Genomics` section. Click on the track name and set each view to `hide` except for `Basewise_Conservation_(phyloP)` which should be set to `full`. Then click `submit`.

*1. Can you see any relationship between these tracks?*


#### **Step 3.3** Use the Table Browser to find correlation between tracks

Sometimes it is hard to tell whether a relationship is real or just perceived. To test how these tracks compare to each other, under the `Tools` menu select `Table Browser`. Select the conservation track you just turned on. Then select the region `chr1:20000000-50000000`. Now click the button `Correlate` and set the window to 20000 bases and click `calculate`.

*1. How do these data relate to each other? Is it what you expected comparing them by eye?*
*2. How would you explain this relationship between SNP density and conservation score?*


### Exercise 4: Perform regression modeling on SNP density data

Given the trends observed in the plot from step 2.2, let's highlight these and get some significance values to go with them.


#### **Step 4.1** Calculate linear model and plot

Add fitted lines using the `geom_smooth` with the method `lm`. Then run a linear regression, grouping the data by feature and save the report.

*1. Are any of the linear models significant at alpha=0.05%? If so, which ones?*

<br><br>

## Submission

Part 1 (2.5 pts):
  * Shell script for running the bedtools commands to create the feature bed files. - 2.5 pts
Part 2. (.57 pts)
  * Shell script for running calculting enrichments - 4.5 pts
    * Using for loops to compare each combination of MAF and feature files using bedtools coverage. - 2 pts
    * Counting number of SNPs and number of bases using awk and storing them in variables - 1 pt
    * Calculate normalized SNP enrichments - 1 pt
    * Writing enrichments along with MAF and feature labels to a file - 0.5 pts
  * Text file with SNP enrichments 0.5 pts 
  * Plot from step 2.4 - 1.5 pts
  * Answer to questions in step 2.4 in README.md - 1 pt

**Total Points: 10**

<br><br>
