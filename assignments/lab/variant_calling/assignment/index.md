# Variant Calling

## Assignment Overview

Today you will perform alignment and variant calling across multiple haploid yeast strains. These strains are the progeny of a cross between a lab strain and a wine strain. The data come from this paper: "[Finding the sources of missing heritability in a yeast cross](http://www.nature.com/nature/journal/v494/n7436/full/nature11867.html)". After performing variant calling, you'll be doing some basic exploratory analysis of the genetic variation you discovered.<br><br>

## Data

##### **Sequencing reads**

As part of this assignment, you'll be working with ten sets of single-end Illumina sequencing reads, each for a different yeast strain. These are in a tarball on Dropbox. You can download the data into your current directory with the following commands:

```
wget https://www.dropbox.com/s/ogx7nbhlhrp3ul6/BYxRM.tar.gz
tar -xvzf BYxRM.tar.gz
```

##### **Reference genome**

You will be aligning reads from your yeast samples to the *Saccharomyces cerevisiae* reference genome. This reference is called _sacCer3_ by the UCSC genome browser, but its name in the NCBI Assembly archive is [R64-1-1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/). **You'll need this info later for snpEff**.

We have provided a single whole genome reference file for you:
```
/Users/cmdb/data/genomes/sacCer3.fa
```

Make a working copy of this file within your weekly homework directory.

As always, if your `.gitignore` is not already set up to ignore these files (the reads and the reference), you should update it so that they are ignored. **You should NOT be uploading any of these files.**<br><br>

## Excercises

There are three exercises in today's assignment:
1. Alignment of sequencing reads to the reference
2. Variant calling and annotation
3. Analysis of annotated variant calls in Python

The first two exercises will occur entirely in the command line. As such, you will be writing and submitting an `alignment-variant_calling.sh` bash script that will run all of the code for exercises 1 and 2.

<span style="color:red;font-weight:bold">NOTE:</span> Some of the steps can take a while. When you've finished running each step, you can comment out that section of your bash script so that it doesn't run every time you run your script. **However**, please uncomment these lines before submitting your script.

Create an empty `alignment-variant_calling.sh` script now.<br><br>

### Exercise 1: Read alignment

During the first exercise, you'll be aligning the reads from each of the 10 yeast strains to the _sacCer3_ reference.<br><br>

#### **Step 1.1**: Index the sacCer3 genome

You'll be using a tool called `bwa` ([BWA manual](http://bio-bwa.sourceforge.net/bwa.shtml)) to perform alignment. Before you can align your sequencing reads, `bwa` needs you to index the sacCer3 genome. Without getting into the the nitty gritty, this essentially means creating a table of contents (and index) that `bwa` can use to quickly find matches between your reads and the reference genome. 

Using `bwa index`, create an index for the `sacCer3.fa` reference.<br><br>

#### **Step 1.2**: Align your reads to the reference

Now that you've indexed the reference, you can align your reads to the reference using `bwa mem`. Because you'll want to run this step the same way on all 10 strains, it makes sense to do this step in a bash `for` loop. Consider this resource for how to write a for loop in a bash script: [bash for loops walkthrough](https://linuxhint.com/bash-for-loop-examples/).

Create a bash `for` loop that loops through each of the 10 samples. For each sample, use `bwa mem` to align the reads to the reference. 

**IT IS VERY IMPORTANT** that you assign each sample a read group during this process, so that individual samples can be distinguished later in Step 2.1. You can do this with the (somewhat cryptic) `-R` flag, which you use to add a line to the header of each output alignment file. An example of a header line you can add with the `-R` flag is `"@RG\tID:Sample1\tSM:Sample1"`. You can replace "Sample1" here with the appropriate sample name for each of your yeast strains.

Perhaps consider the `-t` and `-o` flags as well.<br><br>

#### **Step 1.3**: Format and index your alignments

Now that you've aligned your reads to the reference, you should have 10 `.sam` files, one for each sample. These files contain all of the alignments for each yeast strain. You can see how they're organized with `less -S`, and you can read more about the SAM format [here](https://samtools.github.io/hts-specs/SAMv1.pdf).

These files contain all of the information you need for variant calling, but before you can do that, they'll need to be sorted and indexed (similar to how you indexed the reference in Step 1.1. For both of these tasks you can use the `samtools` program (manual [here](http://www.htslib.org/doc/samtools.html), or you can just run `samtools help`).

**First**, sort each of your `.sam` files using `samtools sort`. You can do this in a new `for` loop in your bash script or, even better, in the same `for` loop you used for alignment. You'll want to output these sorted files as `.bam` files, which contain the same information as the `.sam` file but are compressed.

Perhaps consider the `-O` and `-o` flags when running `samtools sort`. 

**Next**, create an index for each of the resulting sorted `.bam` files using `samtools index`. As before, you can do this in a new `for` loop in your bash script or in the same `for` loop as the previous two steps.

At the end of this step, you should have 10 sorted `.bam` files and their corresponding `.bam.bai` indices.<br><br>

### Exercise 2: Variant calling and annotation

Now that you've aligned the sequencing reads to the reference, you can call genetic variants across the yeast strains.<br><br>

#### **Step 2.1**: Call variants

For variant calling, you'll be using a tool called `freebayes` (manual [here](https://github.com/freebayes/freebayes#usage), or you can just run `freebayes --help`).

Use `freebayes` to identify genetic variants in all of your yeast strains **concurrently** (i.e. you should only be running `freebayes` once will all samples, not for each sample separately). It will output results in Variant Call Format (`.vcf`). You can read more about the VCF format [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

You should consider using the `-f`, `--genotype-qualities`, and `-p` flags. You might like the `-L` flag as well.

**NOTE**: For TAs, running this step took nearly 15 minutes. We expect this step will take a similar amount of time for you, and your computer might make a lot of noise.<br><br>

#### **Step 2.2**: Filter variants based on site quality

Sometimes, `freebayes` will call variants that may not be "real". Luckily, variant callers generally report, for each variant, a "site quality", that describes the probability of that site truly being polymorphic. You can use that quality score to filter out low quality (low confidence) variants.

Note that the "site quality" is not the same as the "genotype quality", which describes the probability that a single sample's genotype at some variant is correct. Both qualities, however, are generally reported on the "Phred" scale (more [here](https://en.wikipedia.org/wiki/Phred_quality_score)).

You can filter out low quality variants using the `vcffilter` tool (documentation [here](https://github.com/vcflib/vcflib/blob/master/doc/vcffilter.md)).

Filter your VCF using `vcffilter` so that you only keep variants whose estimated probability of being polymorphic is greater than 0.99. You should consider how to do this with the `-f` flag. Output your filtered variant calls to a new VCF file.<br><br>

#### **Step 2.3**: Decompose complex haplotypes

Sometimes, especially in regions with complex alignments between samples, you can run into cases where there are multiple possible alternative alleles at the same position. There's nothing inherently "wrong" with this, but it does often complicate downstream analyses.

Luckily, you can use the `vcfallelicprimitives` tool (documentation [here](https://github.com/vcflib/vcflib/blob/master/doc/vcfallelicprimitives.md)) to decompose these more complex haplotypes into more manageable biallelic variants.

Use `vcfallelicprimitives` to decompose complex haplotypes in your filtered VCF from Step 2.2. We suggest using the `-k` and `-g` flags to keep annotations for the variant sites and sample genotypes in your VCF. Output to a new VCF file.<br><br>

#### **Step 2.4**: Annotate variants

Now that you've got these high-quality and nicely behaving variant calls, you want to know what impact these variants might have. Obviously, this is a huge question, and is the basis of all of genetics, but we can get a basic idea of their functional impact on nearby genes (e.g. are they missense, nonsense, etc.) using the `snpEff` tool (the online documentation isn't great, but running `snpEff ann --help` should provide some useful information).

If it wasn't obvious, `snpEff` requires prior annotations (e.g. gene annotations) to work. Have `snpeff` download its database of *Saccharomyces cerevisiae* annotations using the following command (we told you the NCBI ID would be relevant):

```
snpEff download R64-1-1.105
```

Finally, use `snpEff ann` to annotate your VCF with the predicted functional effects that these genetic variants may have. Output to a new (and final) VCF.

For submission purposes, use `head` to grab just the first 100 lines of your final VCF and store this in a new VCF. You will submit this "sample" VCF along with the rest of your assignment. **YOU SHOULD NOT SUBMIT ANY OTHER VCFS; THEY ARE TOO BIG**. Depending on how your `.gitignore` is set up, you may need to do `git add --force <yoursamplevcf.vcf>`.<br><br>

### Exercise 3: Exploratory data analysis

Now that you've discovered variants in these strains and annotated their predicted functional effects, you'd like to do some basic exploratory analysis of the patterns you observe in the VCF. You'll be using Python to create a _single_ nicely formatted and labeled multi-panel plot (e.g., use `subplots` with multiple rows and columns) that explores the following features of the data:
1. The distribution of read depth across sample genotypes
2. The distribution of genotyping quality across samples genotypes
3. The allele frequency spectrum of the discovered variants
4. A summary of the predicted effects of the discovered variants

Each feature will be a single panel in your multi-panel plot.

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

Make sure you label the panel appropriately.<br><br>

#### **Step 3.2**: Genotype quality distribution

Plot a histogram showing the distribution of genotyping quality at each variant across all samples (as before, if you had 10 variants and 5 samples, you'd have 50 data points).

This information can be found in the sample specific FORMAT fields and the end of each line. Check the file header to decide which ID is appropriate. Remember, "genotype quality" is **NOT** the same as "site quality" and is a different part of the VCF line.

Make sure you label the panel appropriately.<br><br>

#### **Step 3.3**: Allele frequency spectrum

Plot a histogram showing the allele frequency spectrum (distribution) of the variants in the VCF (this is a per-variant metric, so with 10 variants and 5 samples, you'd only have 10 data points).

This information is pre-calculated for you and can be found in the variant specific INFO field. Check the file header to decide which ID is appropriate.

Make sure you label the panel appropriately.<br><br>

#### **Step 3.4**: Predicted effects

Create a barplot showing the predicted effect(s) of the variants in the VCF (i.e. a bar for each "type" of effect showing the number of variants with that effect).

In Python, produce a nicely formatted and labeled multi-panel plot (e.g., use `subplots` with multiple rows and columns) describing your variants.<br /><br />Explore each of the following characteristics of the variant genotypes called across all ten yeast samples. (Each characteristic will be a subplot in the multi-panel plot).

This information was added to the VCF by `snpEff` and can be found in the variant specific INFO field. Check the file header to decide which ID is appropriate and how to parse the information.

**NOTE**: We encourage you to consider every possible effect for each variant, but feel free to just grab the first one.

Make sure you label the panel appropriately.<br><br>

## Submission

1. Bash script that performs read alignment and variant calling/filtering/annotation (**3 points**)
 * Read alignment (**1 point**)
 * Variant calling (**1 point**)
 * Filtering and annotation (**1 point**)
2. VCF file containing the first 100 lines of your filtered/annotated VCF from Step 2.4 (**1 point**)
 * **DO NOT SUBMIT ANY OF THE OTHER VCFS; THEY ARE TOO BIG**
3. Python script that runs exploratory analysis of VCF from Step 2.4 (**4 points**)
   * Code to parse VCF file (**2 points**)
   * Code to produce multi-panel plot (**2 points**)
4. Nicely formatted and labeled multi-panel plot showing summaries of exploratory analysis (**2 points**)

**Total Points: 10**

<br><br>
