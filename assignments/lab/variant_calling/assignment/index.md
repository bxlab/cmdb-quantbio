# Assignment 3: Variant Calling
Assignment Date: Friday, Sept. 23, 2022 <br>
Due Date: Friday, Sept. 30, 2022 @ 1:00pm ET <br>

**Please do a `git pull` within your `~/cmdb-quantbio/` directory**

## Assignment Overview

Today we will perform <i>de novo</i> identification of variants in multiple haploid yeast strains. These strains are the progeny of a cross between a lab strain and a wine strain. The data come from [Finding the sources of missing heritability in a yeast cross](http://www.nature.com/nature/journal/v494/n7436/full/nature11867.html).

## Data

##### **Sequencing reads**

The following zip file contains ten sets of single-end Illumina sequencing reads, each for a different yeast strain.

```
/Users/cmdb/cmdb-quantbio/assignments/lab/variant_calling/extra_data/BYxRM.tar.gz
```

Make a working copy of this file within your weekly homework directory and then unzip the file using the following command.

```
tar -xvzf BYxRM.tar.gz
```

##### **Reference genome**

You will be aligning reads from your yeast samples to the *Saccharomyces cerevisiae* reference genome. This reference is called _sacCer3_ by the UCSC genome browser, but its name in the NCBI Assembly archive is [R64-1-1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/). **You need this info later for snpEff**.

We have provided a single whole genome reference file for you. (all of the chromosomes are in the FASTA same file as separate entries)

```
/Users/cmdb/data/genomes/sacCer3.fa
```

Make a working copy of this file within your weekly homework directory. If you try to work with the file at this original path, you will encounter directory write permission errors.

## Goal of Assignment

Your task now is to identify genetic variants in these samples and annotate them with their predicted functional effects.

---

### Submitting your assignment

If you do not write a bash script, keep track of all your command-line code in a `README.md` file.

Push

* all scripts (Python and/or command line),
* your `README.md` file (if applicable),
* a nicely formatted multi-panel plot,
* and **ONLY the first 1000 lines** of your filtered, decomposed, and annotated VCF to your `answers` Github repository. **Do not push any other raw data to Github, and do not push the full VCF!**
  * You will need to do `git add --force <yourvcfilename.vcf>`


### Bash scripts

Since several parts (Steps 2, 3a, and 3b) of this assignment will involve running the same code on each of your 10 yeast samples, you may want to write a bash script that will automate this process for you. However, this is *not* required.

(Note: The sample names are the names of provided sequencing read fastq files)

Consider this resource for how to write a for loop in a bash script: [bash for loops walkthrough](https://linuxhint.com/bash-for-loop-examples/)

### Step 1: Index the sacCer3 genome with `bwa index`

Before you can align your sequencing reads, you first need to index the sacCer3 genome using `bwa index`.

### Step 2: Alignment with `bwa mem`

Align your reads against the sacCer3 reference genome using `bwa mem`. You will want to run this command 10 separate times (one for each fastq file).

**IT IS VERY IMPORTANT** that you assign each sample a read group during this process, so that individual samples can be distinguished in Step 4. You can do this with the (somewhat cryptic) `-R` flag, which you use to add a line to the header of each output alignment file. An example of a header line you can add with the `-R` flag is `"@RG\tID:Sample1\tSM:Sample1"`. You can replace "Sample1" here with the appropriate sample name for each fastq file.

Perhaps consider the `-t` and `-o` flags as well.

[BWA manual](http://bio-bwa.sourceforge.net/bwa.shtml)

### Step 3a: Create a sorted bam file with `samtools sort`, for input to variant callers

You will need to create a sorted bam file using `samtools sort`. Again, you will want to run this command 10 separate times (one for each sample output from `bwa mem`)

Perhaps consider the `-O` and `-o` flags.

[samtools manual](http://www.htslib.org/doc/samtools-sort.html)


### Step 3b: Create an index for each sorted bam file with `samtools index`, for input to variant callers

You will want to create an index file for each sorted bam file using `samtools index`. Consider the `-b` flag and name the output files with the `.bam.bai` extension to signify that they are indexed bam files. Again, you will want to run this command 10 separate times (one for each sample output from `samtools sort`)

### Step 4: Variant calling with `freebayes`

Use `freebayes` to identify genetic variants in all of your yeast strains concurrently. It will output results in Variant Call Format (VCF). You should consider using the `-f`, `--genotype-qualities`, and `-p` flags. You might like the `-L` flag as well.

[freebayes documentation](https://github.com/ekg/freebayes#getting-the-best-results)

Note: For TAs, running this step took nearly 15 minutes. This step will take a similar amount of time for you, and your computer might make a lot of noise.

### Step 5: Filter variants based on genotype quality using `vcffilter`

Filter your VCF using `vcffilter` so that you only keep variants whose estimated probability of being polymorphic is greater than 0.99. You should consider how to do this with the `-f` flag. The `freebayes` documentation will be helpful here, as well as [this vcffilter info](https://github.com/vcflib/vcflib#vcffilter).

### Step 6: Decompose complex haplotypes using `vcfallelicprimitives`

Use `vcfallelicprimitives` to decompose complex haplotypes within the filtered variants file. We suggest using the `-k` and `-g` flags to keep annotations for the variant sites and sample genotypes in your VCF.

You can reference [vcfallelicprimitives documentation 1](https://github.com/vcflib/vcflib#vcfallelicprimitives) and [vcfallelicprimitives documentation 2](https://janis.readthedocs.io/en/latest/tools/bioinformatics/ekg/vcfallelicprimitives.html).

### Step 7: Variant effect prediction with `snpeff ann`

First, you'll need to downgrade/re-install `snpeff`:

```
conda install snpeff=5.0 -y
```

Then fetch the appropriate yeast reference database:

```
snpeff download R64-1-1.99
```

Then, use `snpeff ann` to annotate your VCF with the predicted functional effects that these genetic variants may have.

We recommend *not* Googling the `snpeff` documenation. It will tell you to use `java -jar snpEff.jar`, which you should not. The help option for `snpeff ann`'s command-line tool is 100 times better.

### Step 8: Exploratory data analysis through plotting

In Python, produce a nicely formatted and labeled multi-panel plot (e.g., use `subplots` with multiple rows and columns) describing your variants.<br /><br />Explore each of the following characteristics of the variant genotypes called across all ten yeast samples. (Each characteristic will be a subplot in the multi-panel plot).

  * The read depth distribution of variant genotypes (histogram)
      * This information can be found in the sample specific FORMAT field for each variant/line. Check the file header to decide which ID is appropriate.
  * The quality distribution of variant genotypes (histogram)
      * This information can be found in the sample specific FORMAT field for each variant/line. Check the file header to decide which ID is appropriate.
  * The allele frequency spectrum of your identified variants (histogram)
      * This information is pre-calculated for you and can be found in the variant specific INFO field. Check the file header to decide which ID is appropriate.
  * A summary of the predicted effect(s) of each variant as determined by snpEff (barplot)
      * This information was added to the VCF by snpEff and can be found in the variant specific INFO field. Check the file header to decide which ID is appropriate and how to parse the information.
      * We encourage you to consider every possible effect for each variant, but feel free to just grab the first one.

You may find it helpful to reference [this page](https://pcingola.github.io/SnpEff/se_inputoutput/) of the `snpeff` manual, which describes the format of its output VCF.

### Submit!

Push all scripts, a record of your command line commands (if applicable), your multi-panel plot, and **ONLY the first 1000 lines** of your filtered, decomposed, and annotated VCF to your `qbb2022-answers` repo. **Do not push any other raw data to Github, and do not push the full VCF!**
