# Variant Calling

## Assignment Overview

Modern genome sequencing technologies produce sequencing "reads" which represent subsequences of the fragmented DNA molecules isolated in an experiment. One of the first steps in nearly all genomic analyses is to align (or "map") the reads to a reference genome and identify sites where your samples vary from the reference. This lab is designed to practice and explore this basic approach.

To accomplish this goal, we will be working with Illumina short-read sequencing data generated in a study where a lab strain of a model yeast species (Saccharomyces cerevisiae) was crossed with a wine strain as part of a quantitative genetic experiment to map associations with phenotypes (paper: "[Finding the sources of missing heritability in a yeast cross](http://www.nature.com/nature/journal/v494/n7436/full/nature11867.html)"). The diploid offspring from this cross were sporulated, and one (haploid) spore from each of the resulting tetrads was sequenced. This means that the data you will be working with are from haploid genomes that we expect will be mosaics of the two original strains (i.e., will possess "tracts" of ancestry from one strain versus the other, with relatively few transitions).

## Data

##### **Sequencing reads**

We have randomly selected ten samples from the study and placed them in a compressed file on Dropbox. You can download the data into your current directory and decompress them with the following commands:

```
wget https://www.dropbox.com/s/ogx7nbhlhrp3ul6/BYxRM.tar.gz
tar -xvzf BYxRM.tar.gz
```

##### **Reference genome**

You will be aligning reads from your yeast samples to the *Saccharomyces cerevisiae* reference genome. This reference is called _sacCer3_ by the UCSC genome browser, but its name in the NCBI Assembly archive is [R64-1-1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/).

We have provided a single whole genome reference file for you:
```
/Users/cmdb/data/genomes/sacCer3.fa
```

Make a working copy of this file within your weekly homework directory.

As always, if your `.gitignore` is not already set up to ignore these files (the reads and the reference), you should update it so that they are ignored. **You should NOT be uploading any of these files.**<br><br>


## Excercises

Remember to record all of your work in a shell script! You can answer the questions using comments within your script. Show how you obtained your answers.

### Exercise 1: Get to know your data

Take a glance at the first file, `A01_09.fastq`. These are single-end whole-genome DNA sequencing data from an Illumina sequencer. Use your knowledge of FASTQ files (see Mike Sauria's previous lecture) and the UNIX commands you have learned from previous classes to answer the following questions:

1. How long are the sequencing reads?

77 - 1 (newline) = 76 bp per read

2. How many reads are present within the file?

3. Given your answers to 1 and 2, as well as knowledge of the length of the S. cerevisiae reference genome, what is the expected average depth of coverage?

4. While you do not need to repeat for all samples, looking at the size of the files can give us information about whether we have similar amounts of data from other samples. Use the `du` command to check the file sizes of the rest of the samples. Which sample has the largest file size (and what is that file size, in megabytes)? Which sample has the smallest file size (and what is that file size, in megabytes)?

5. Run the program FastQC on your samples (with default settings). Open the HTML report for sample `A01_09`. What is the median base quality along the read? How does this translate to the probability that a given base is an error? Do you observe much variation in quality with respect to the position in the read? 

### Exercise 2: Map your reads to the reference genome

#### **Step 2.1**: Index the sacCer3 genome

You'll be using a tool called `bwa` ([BWA manual](http://bio-bwa.sourceforge.net/bwa.shtml)) to perform alignment. Before you can align your sequencing reads, `bwa` needs you to index the sacCer3 genome. Without getting into the the nitty gritty, this essentially means creating a table of contents (and index) that `bwa` can use to quickly find matches between your reads and the reference genome. 

Using `bwa index`, create an index for the `sacCer3.fa` reference.<br><br>

#### **Step 2.2**: Align your reads to the reference

Now that you've indexed the reference, you can align your reads to the reference using `bwa mem`. Because you'll want to run this step the same way on all 10 strains, it makes sense to do this step in a bash `for` loop (check your notes from last week!).

Create a bash `for` loop that loops through each of the 10 samples. For each sample, use `bwa mem` to align the reads to the reference. 

**IT IS IMPORTANT** that you assign each sample a read group during this process, so that individual samples can be distinguished later in Step 2.1. You can do this with the (somewhat cryptic) `-R` flag, which you use to add a line to the header of each output alignment file. An example of a header line you can add with the `-R` flag is `"@RG\tID:Sample1\tSM:Sample1"`. You can replace "Sample1" here with the appropriate sample name for each of your yeast strains.

Perhaps consider the `-t` flag as well.<br><br>

#### **Step2.3**: Sanity check your alignments

Now that you've aligned your reads to the reference, you should have 10 `.sam` files, one for each sample. These files contain all of the alignments for each yeast strain. You can see how they're organized with `less -S`, and you can read more about the SAM format [here](https://samtools.github.io/hts-specs/SAMv1.pdf).

Using various Unix commands, answer the following questions about the `A01_09` SAM file:

1. How many chromosomes are in the yeast genome?

2. How many total read alignments are recorded in the SAM file?

3. How many of the alignments are to loci on chromosome III?


#### **Step 2.4**: Format and index your alignments

These files contain all of the information you need for variant calling, but before you can do that, they'll need to be sorted and indexed (similar to how you indexed the reference in Step 1.1. For both of these tasks you can use the `samtools` program (manual [here](http://www.htslib.org/doc/samtools.html), or you can just run `samtools help`).

**First**, sort each of your `.sam` files using `samtools sort`. You can do this in a new `for` loop in your bash script or, even better, in the same `for` loop you used for alignment. You'll want to output these sorted files as `.bam` files, which contain the same information as the `.sam` file but are compressed.

Perhaps consider the `-O` and `-o` flags when running `samtools sort`. 

**Next**, create an index for each of the resulting sorted `.bam` files using `samtools index`. As before, you can do this in a new `for` loop in your bash script or in the same `for` loop as the previous two steps.

At the end of this step, you should have 10 sorted `.bam` files and their corresponding `.bam.bai` indices.<br><br>

#### **Step 2.5**: Visualize your alignments

Open IGV and set the SacCer3 genome as the reference. Load the sample `A01_09` in IGV ("File" -> "Load from File..."). Zoom in far enough to see the reads and scan through some alignments.

1. Does the depth of coverage appear to match that which you estimated in Step 1.3? Why or why not?

2. Set your window to "chrI:113,113-113,343" (paste that string in the search bar and click enter). How many SNPs do you observe in this window? Are there any SNPs about which you are uncertain? Explain your answer.

3. Set your window to "chrIV:825,548-825,931". What is the position of the SNP in this window? Does this SNP fall within a gene?

### Exercise 3: Variant calling and annotation

Now that you've aligned the sequencing reads to the reference genome, you can call genetic variants across the yeast strains. The most widely used program for this purpose is called GATK. I am confident that you could figure out how to use it given enough time, but we would spend the whole class debugging esoteric details of this specific program. To save you the effort, I therefore went ahead and called variants on the BAM files another program called FreeBayes, followed by a bit of quality filtering and adjustment of formatting:

```
ls *.bam > bamListFile.txt

freebayes -f sacCer3.fa -L bamListFile.txt --genotype-qualities -p 1 > unfiltered.vcf

vcffilter -f "QUAL > 20" unfiltered.vcf > filtered.vcf

vcfallelicprimitives -kg filtered.vcf > decomposed_filtered.vcf
```

That file can be obtained here: https://www.dropbox.com/scl/fi/cuk8g4p4wu5atelh0y9y1/biallelic_decomposed_filtered.vcf?rlkey=a3dshq66t0wqycgfzxnqn54gg&dl=0


### Exercise 3: Exploratory data analysis

Now that you've discovered variants in these strains and annotated their predicted functional effects, you'd like to do some basic exploratory analysis of the patterns you observe in the VCF. You will be creating figures that explore the following features of the data:
1. The distribution of read depth across sample genotypes
2. The distribution of genotyping quality across samples genotypes
3. The allele frequency spectrum of the discovered variants

Create an empty `variation_analysis.py` script now where you'll be doing the analyses in the next steps.<br><br>

#### **Step 3.0**: Parse the VCF file

For these analyses, you'll have to parse the filtered/biallelic VCF you generated in Step 2.4. If you can find a python library that handles VCF parsing, you're welcome to use it, but it may be easier to simply use the following structure:

```
for line in open(<vcf_file_name>):
    if line.startswith('#'):
        continue
    fields = line.rstrip('\n').split('\t')

    # grab what you need from `fields`
```

Each of the following analyses needs information from a different field/column from the VCF. As such, you'll be grabbing multiple pieces of information from each line, so if you do use this structure, instead of looping through the VCF separately for each analysis, try to use just a single loop across all analyses.

Remember, you can read more about the VCF file format [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).<br><br>

#### **Step 3.1**: Read depth distribution

Plot a histogram showing the distribution of read depth at each variant across all samples (e.g. if you had 10 variants and 5 samples, you'd have 50 data points).

This information can be found in the sample specific FORMAT fields and the end of each line. Check the file header to decide which ID is appropriate.

Make sure you label the plot appropriately. Use the `scale_x_log10()` option to ggplot2 to scale the x-axis.

Interpret this figure in two or three sentences in your own words. Does it look as expected? Why or why not? 

Bonus: what is the name of this distribution?

<br><br>


#### **Step 3.2**: Allele frequency spectrum

Plot a histogram showing the allele frequency spectrum (distribution) of the variants in the VCF (this is a per-variant metric, so with 10 variants and 5 samples, you'd only have 10 data points). 

This information is pre-calculated for you and can be found in the variant specific INFO field. Check the file header to decide which ID is appropriate.

Make sure you label the panel appropriately. Set `bindwidth=0.025` to avoid binning artifacts.

Interpret this figure in two or three sentences in your own words. Does it look as expected? Why or why not?

Bonus: what is the name of this distribution?

<br><br>


## Submission

1. Bash script that performs the explorations of the FASTQ file for exercise 1 (**2.5 points (0.5 per question)**).

2. Bash script that performs the alignments, formatting, indexing, and answers to questions from exercise 2 (**2.5 points (0.5 point per step)**).

3. Python script to produce the output necessary for plots in Step 3 (**3 points**).

4. R script to take the output from step 3 and generate figures, which should also be uploaded (**2 points**).


**Total Points: 10**

<br><br>
