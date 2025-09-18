# Sequence Alignment

## Overview

Learning Objectives
- Work with [short-read](https://wikipedia.org/wiki/Massive_parallel_sequencing) ([Illumina](https://wikipedia.org/wiki/Illumina_dye_sequencing)) and [long-read](https://wikipedia.org/wiki/Third-generation_sequencing) ([Nanopore](https://wikipedia.org/wiki/Nanopore_sequencing)) sequencing datasets
- Align sequencing data using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2), [minimap2](https://github.com/lh3/minimap2), and [HISAT2](https://daehwankimlab.github.io/hisat2)
- Explore [sequence alignments](https://samtools.github.io/hts-specs) and visualize using [IGV](https://igv.org/doc/desktop)
- Practice constructing Unix commands, creating [Bash scripts](https://swcarpentry.github.io/shell-novice/06-script.html), and summarizing data using Python

## Instructions

- Work in `~/qbXX-answers/weekX`
- Document your Unix commands in a `README.md`
- `git push` after each exercise
- Strive to organize your directory e.g.

    ```
    /Users/cmdb/qb25-answers/week2
    ├── README.md
    ├── genomes
    │   └── sacCer3.fa
    ├── rawdata
    │   └── ERR8562478.fastq
    └── variants
        ├── A01_01.bam
        ├── ...
        └── A01_06.sam
    ```

## Preparation 

Download [BYxRM dataset](http://genomics-pubs.princeton.edu/YeastCross_BYxRM) and unpack using `tar`

```
cd ~/Data
wget -O BYxRM.tar.gz "https://www.dropbox.com/scl/fi/ggmmx6ceqni3e2wnljn0t/BYxRM.tar.gz?rlkey=rfkthvs49pa56y9eo3ixyq9ep&dl=0"
tar xf BYxRM.tar.gz
```

Confirm that `ls -l BYxRM` matches

```
-rw-r--r--   1 cmdb  staff  23729508 May  6  2024 BYxRM_GenoData.txt
-rw-r--r--   1 cmdb  staff    797659 May  6  2024 BYxRM_PhenoData.txt
drwxr-xr-x  98 cmdb  staff      3136 Jul 30 15:13 fastq
```

Skim genotype calls in BYxRM_GenoData.txt using `less`

```
marker                 A01_01  A01_02
27915_chr01_27915_T_C  R       B
28323_chr01_28323_G_A  R       B
28652_chr01_28652_G_T  R       B
29667_chr01_29667_C_A  R       B
```

## Exercises

1. Align short sequencing reads using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner)

    [Build](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) reference genome index (only one time per reference)

    ```
    cd genomes
    cp ~/Data/References/sacCer3/sacCer3.fa.gz .
    gunzip sacCer3.fa.gz
    bowtie2-build sacCer3.fa sacCer3
    ```

    [Map](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example) reads to genome and [sort/index](https://en.wikipedia.org/wiki/SAMtools#Examples) with samtools (as many times as samples)

    ```
    cd variants
    bowtie2 -p 4 -x ___/sacCer3 -U ___/A01_01.fq.gz > A01_01.___
    samtools sort -o ___ ___
    samtools index ___
    samtools idxstats ___ > A01_01.idxstats
    ```

    [Visualize](https://igv.org/doc/desktop/#UserGuide/navigation) reads using IGV.app

    - Switch to "S. cerevisiae (sacCer3)" genome in the top left menu
        - If not available, in the Genomes menu select "Download Hosted Genome..."
    - Load A01_01.bam
    - View region chrI:27,000-32,000
    - View genotype calls in BYxRM_GenoData.txt
        - chrI:27915
        - chrI:28323
        - chrI:28652
        - chrI:29667
    - File > Save PNG Image ... of chrI:27,000-32,000 as A01_01.png

    Submit the following

    - README.md -- One `bowtie2` and three `samtools` commands
    - A01_01.idxstats -- Output of `samtools idxstats`
    - A01_01.png -- Image of chrI:27,000-32,000 region

1. Build a workflow to process the first six samples using a [Bash script](https://swcarpentry.github.io/shell-novice/06-script.html)

    Create `map-reads.sh` by completing this template

    ```
    #!/bin/bash

    for sample in # list prefixes e.g. A01_01 A01_02
    do
        echo "***" $sample
        # mapping command e.g. bowtie2 path/to/$sample.fq.gz
        # sort command
        # index command
    done
    ```

    Visualize haplotype using IGV.app

    - Load all six .bam files
    - View region chrI:27,000-32,000
    - Click - to zoom out twice
    - Save PNG Image ... as A01_01-06.png
    - Describe in README.md how this visualization compares to [haplotypes](https://www.genome.gov/genetics-glossary/haplotype) in BYxRM_GenoData.txt

    Submit the following

    - README.md -- Your description of the visualization
    - map-reads.sh -- Bash script to process six samples along with comment on visualization
    - A01_01-06.png -- Image of chrI:27,000-32,000 region

1. Summarize sequence alignments using Python

    Create `summarize-sam.py` to count alignments to each chromosome

    - Process line by line a [.sam file](https://samtools.github.io/hts-specs/SAMv1.pdf) specified as the first command line argument
    - Skip header lines that begin with `@`
    - Use a dictionary and count how many alignments there are for each chromosome (and unmapped) as reported in the `RNAME` field
    - Print each dictionary key, value pair in the default order returned by .keys()
    - Comfirm that you get the same results as `samtools idxstats`
        - If you want your output ordered similarly, sort the input by `samtools sort -o reads.sorted.sam reads.sam`

    Extend summarize-sam.py to examine mismatches per alignment

    - Count how many times each `NM:i:count` [SAM tag](https://samtools.github.io/hts-specs/SAMtags.pdf) occurs
    - Note that `NM` is not always in the same column so use a for loop to go through the fields after splitting a line
    - For the dictionary key, remove the `NM:i:` by slicing and convert to `int()`
    - Print the dictionary keys in numerical order by first using the `sorted()` function e.g.

    Submit the following

    - summarize-sam.py -- Python script to count chromosomes and mismatches
    - summarize-sam.txt -- Output of summarize-sam.py on A01_01.sam

1. Align long sequencing reads using [minimap2](https://github.com/lh3/minimap2)

    Find yeast Nanopore dataset at [SRA](https://www.ncbi.nlm.nih.gov/sra)

    - Read the project description at https://www.ncbi.nlm.nih.gov/bioproject/PRJEB50706 
    - Scroll down and click on the 288 [SRA Experiments](https://www.ncbi.nlm.nih.gov/sra/docs/submitmeta/#sra-metadata-experiment)
    - Use the left hand column to filter for Oxford Nanopore
    - Click on ERX8178858 (should be #10 on the list)
    - Click on the Run number (ERR8562478) to see the dataset size, average read length, taxonomy analysis

    Fetch reads using [sratoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/#download-sequence-data-files-usi)

    - Change into your `week2/rawdata` directory
    - Download dataset using `fasterq-dump -p ERR8562478`
        - If not available, run `conda install sra-tools`
    - Confirm that `ERR8562478.fastq` has 404400 lines

    Map reads to genome

    - Create and work in `week2/longreads`
    - Run `minimap2` with the following [arguments](https://github.com/lh3/minimap2#general-usage)
        - If not available, run `conda install minimap2`
        - Input/Output option to output in the SAM format
        - Preset option to map Nanopore reads
        - sacCer3.fa as the target.fa
        - ERR8562478.fastq as the query.fa
        - Redirecting output to a file named `longreads.sam`
    - Sort .sam into .bam file and index
    - Save the output of samtools idxstats as longreads.idxstats

    Visualize using IGV.app

    - Open in IGV and visualize chrIII:148,000-153,000
    - Click the View menu, select Preferences..., switch to the [Third Gen](https://www.pacb.com/blog/igv-3-improves-support-pacbio-long-reads/) tab, turn on "Hide indels < show indel threshold" and set to "100" 
    - Note the structural variants ([transposon insertions](https://pubmed.gov/17173485)) in these samples relative to the reference genome indicated by purple rectangles ~5800 nt long
    - Save PNG Image ... as chrIII-150kb.png

    Submit the following

    - README.md -- One `minimap2` command
    - longreads.idxstats -- Output of `samtools idxstats`
    - chrIII-150kb.png -- Image of chrIII:148,000-153,000 region

1. Align RNA-seq reads using [HISAT2](https://daehwankimlab.github.io/hisat2)

    - Download SRR10143769 and confirm that it has 11670744 lines
    - Build reference genome index using `hisat2-build sacCer3.fa sacCer3`
    - Create and work in `week2/rna`
    - Map reads to genome using `hisat2` (command similar to `bowtie2`)
    - Sort, index, and open in IGV
    - Save an image of a region with three active genes that have >20 reads per gene
    - Describe in README.md what part of the genes appear to have the most coverage

    Submit the following

    - README.md -- One `hisat2` command and your description
    - rna.png -- Image showing three active genes with RNA expression

## Grading

- 1 pt -- Align short sequencing reads
- 3 pt -- Build a workflow
- 3 pt -- Summarize sequence aligments
- 2 pt -- Align long sequencing reads
- 1 pt -- Align RNA-seq reads 
