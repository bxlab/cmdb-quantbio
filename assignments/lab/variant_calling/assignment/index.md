# Assignment 3: Variant Calling
Assignment Date: Friday, Sept. 23, 2022 <br>
Due Date: Friday, Sept. 30, 2022 @ 1:00pm ET <br>


## Lecture

**Slides** are available here: `~/cmdb-quantbio/assignments/lab/variant_calling/slides_asynchronous_or_livecoding_resources/20210914_qblab_variant_calling.pptx`


## Assignment Overview

Today we will perform <i>de novo</i> identification of variants in multiple haploid yeast strains. These strains are the progeny of a cross between a lab strain and a wine strain. The data come from [Finding the sources of missing heritability in a yeast cross](http://www.nature.com/nature/journal/v494/n7436/full/nature11867.html).

## Data

##### **Sequencing reads**

The following zip file contains ten sets of single-end Illumina sequencing reads, each for a different yeast strain.

```
wget "https://github.com/bxlab/qbb2021/raw/main/week2/BYxRM_subset.tar.xv"
tar -xvzf BYxRM_subset.tar.xv
```

##### **Reference genome**

You will be aligning reads from your yeast samples to the *Saccharomyces cerevisiae* reference genome. This reference is called _sacCer3_ by the UCSC genome browser, but its name in the NCBI Assembly archive is [R64-1-1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/). **You need this info later for snpEff**.

When you unzip the reference genome file from UCSC, it will have separate fasta files for each chromosome. You should combine these files into a single whole-genome reference by using the code below:

```
wget "http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz"
tar -xvzf chromFa.tar.gz
cat chr*.fa > sacCer3.fa
rm chr*.fa
```

Your task now is to identify genetic variants in these samples and annotate them with their predicted functional effects.

---

### Submitting your assignment

If you do not write a bash script, keep track of all your command-line code in a `txt` or `md` file.

Push all scripts (Python and/or command line), your `txt` or `md` file (if applicable), a nicely formatted multi-panel plot, and **ONLY the first 1000 lines** of your filtered, decomposed, and annotated VCF to your `answers` Github repository. **Do not push any other raw data to Github, and do not push the full VCF!**



### Bash scripts

Since several parts of this assignment will involve running the same code on each of your 10 yeast samples, you may want to write a bash script that will automate this process for you. However, this is *not* required.

### Step 1: Index the sacCer3 genome with `bwa index`

Before you can align your sequencing reads, you first need to index the sacCer3 genome.

### Step 2: Alignment with `bwa mem`

Align your reads against the sacCer3 reference genome.

**IT IS VERY IMPORTANT** that you assign each sample a read group during this process, so that individual samples can be distinguished in Step 4. You can do this with the (somewhat cryptic) `-R` flag, which you use to add a line to the header of each output alignment file. An example of a header line you can add with the `-R` flag is `"@RG\tID:Sample1\tSM:Sample1"`.

Perhaps consider the `-t` and `-o` flags as well.

[BWA manual](http://bio-bwa.sourceforge.net/bwa.shtml)

### Step 3: Create a sorted bam file with `samtools`, for input to variant callers

Perhaps consider the `-O` and `-o` flags.

[samtools manual](http://www.htslib.org/doc/samtools-sort.html)

### Step 4: Variant calling with `freebayes`

Use `freebayes` to identify genetic variants in all of your yeast strains concurrently. It will output results in Variant Call Format (VCF). You should consider using the `-f`, `--genotype-qualities`, and `-p` flags.

[freebayes documentation](https://github.com/ekg/freebayes#getting-the-best-results)

Note: This step will take a few minutes, and your computer might make a lot of noise.

### Step 5: Filter variants based on genotype quality using `vcffilter`

Filter your VCF so that you only keep variants whose estimated probability of being polymorphic is greater than 0.99. You should consider how to do this with the `-f` flag. The `freebayes` documentation will be helpful here, as well as [this vcffilter info](https://github.com/vcflib/vcflib#vcffilter).

### Step 6: Decompose complex haplotypes using `vcfallelicprimitives`

We suggest using the `-k` and `-g` flags to keep annotations for the variant sites and sample genotypes in your VCF.

You can reference [vcfallelicprimitives documentation 1](https://github.com/vcflib/vcflib#vcfallelicprimitives) and [vcfallelicprimitives documentation 2](https://janis.readthedocs.io/en/latest/tools/bioinformatics/ekg/vcfallelicprimitives.html).

### Step 7: Variant effect prediction with `snpeff ann`

First, fetch the appropriate yeast reference database:

```
snpeff download R64-1-1.99
```

Then, use `snpeff ann` to annotate your VCF with the predicted functional effects that these genetic variants may have.

We recommend *not* Googling the `snpeff` documenation. It will tell you to use `java -jar snpEff.jar`, which you should not. The help option for `snpeff ann`'s command-line tool is 100 times better.

### Step 8: Exploratory data analysis through plotting

In Python, produce a nicely formatted and labeled multi-panel plot describing your variants.<br /><br />Explore each of the following characteristics of the variant genotypes called across all ten yeast samples. (Each characteristic will be a subplot in the multi-panel plot).

  * The read depth distribution of variant genotypes (histogram)
      * This information can be found in the sample specific FORMAT field for each variant/line. Check the file header to decide which ID is appropriate.
  * The quality distribution of variant genotypes (histogram)
      * This information can be found in the sample specific FORMAT field for each variant/line. Check the file header to decide which ID is appropriate.
  * The allele frequency spectrum of your identified variants (histogram)
      * This information is pre-calculated for you and can be found in the variant specific INFO field. Check the file header to decide which ID is appropriate.
  * A summary of the predicted effect(s) of each variant as determed by snpEff (barplot)
      * This information was added to the VCF by snpEff and can be found in the variant specific INFO field. Check the file header to decide which ID is appropriate and how to parse the information.
      * We encourage you to consider every possible effect for each variant, but feel free to just grab the first one.

You may find it helpful to reference [this page](https://pcingola.github.io/SnpEff/se_inputoutput/) of the `snpeff` manual, which describes the format of its output VCF.

### Submit!

Push all scripts, a record of your command line commands (if applicable), your multi-panel plot, and **ONLY the first 1000 lines** of your filtered, decomposed, and annotated VCF to your `qbb2021-answers` repo. **Do not push any other raw data to Github, and do not push the full VCF!**

---

"Got anything else?" you ask. Of course we do. If you have time, try your hand at the advanced exercises:

## Advanced Exercises: Coverage Simulation & de Bruijn Graphs

### Submitting your assignment

Push all scripts, a `txt` or `md` file with your answers, and plots to your `qbb2021-answers` Github repository.

### Question 1. Coverage simulator

- Q1a. How many 100bp reads are needed to sequence a 1Mbp genome to 5x coverage?

- Q1b. Using Python, simulate sequencing 5x coverage of a 1Mbp genome and plot a histogram of the coverage. Overlay the histogram with a Poisson distribution with lambda = 5.
    - Note you do not need to actually output the sequences of reads - you can randomly sample positions in the genome and continually record the "coverage" you get from this sampling.
    - You do not need to consider the strand of reads.
    - Assume that sequencing reads have a uniform random probability of starting at each possible position (1 through 999,900).
    - You can record the coverage in an array of 1M positions.

- Q1c. Using the histogram from 1b, how much of the genome has not been sequenced (has 0x coverage)? How well does this match Poisson expectations?

- Q1d. Now repeat the analysis with 15x coverage:
    1. Simulate the appropriate number of reads
    2. Make a histogram
    3. Overlay a Poisson distribution with lambda = 15
    4. Compute the number of bases with 0x coverage
    5. Evaluate how well it matches the Poisson expectation

### Question 2. de Bruijn graph construction

- Q2a. Draw, with Python, the de Bruijn graph for the following reads using k=3. Assume all reads are from the forward strand, no sequencing errors, and complete coverage of the genome.

```
ATTCA
ATTGA
CATTG
CTTAT
GATTG
TATTT
TCATT
TCTTA
TGATT
TTATT
TTCAT
TTCTT
TTGAT
```

- Q2b. Using the k-mers from Q2a, ssume that the maximum number of occurrences of any 3-mer in the actual genome is 4. Write out one possible genome sequence.

- Q2c. What would it take to fully resolve the genome?
