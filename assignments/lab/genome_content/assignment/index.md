# Genome Content

## Overview

Learning Objectives
- Explore [genome features](https://en.wikipedia.org/wiki/Human_genome) such as genes, chromatin, and variants
- Analyze data using [UCSC Genome Browser](https://genome.ucsc.edu) and [BEDTools](https://bedtools.readthedocs.io)
- Practice running programs in Unix and plotting data in R

## Instructions

Document your answers in `~/qbXX-answers/weekX`

`git push` after each exercise and **do not wait** until the end of the lab

## Exercises

1. Plot gene density across each chromosome

    Create hg19-main.chrom.sizes with information about the 25 main chromosomes

    - Use `wget` to obtain `hg19.chrom.sizes` from https://hgdownload.soe.ucsc.edu/downloads.html under "Genome sequence files"
    - View the file using `less` to see primary and secondary contigs
    - Use the following command to exclude haplotypes, unplaced, and random contigs and rename chrM to chrMT
        ```
        grep -v _ hg19.chrom.sizes | sed 's/M/MT/' > hg19-main.chrom.sizes
        ```

    Create hg19-1mb.bed with 1 mb intervals across the hg19 assembly

    - Type `bedtools makewindows` to view syntax and documentation
        ```
        Usage: bedtools makewindows [OPTIONS] [-g <genome> OR -b <bed>]
        [ -w <window_size> OR -n <number of windows> ]
        ```
    - Run `bedtools makewindows` and save the output using `>`
    - Confirm that hg19-1mb.bed has 3,114 lines

    Create hg19-kc.bed with one transcript per gene (aka [knownCanonical](https://genome.ucsc.edu/FAQ/FAQgenes.html#singledownload))

    - Navigate to https://genome.ucsc.edu and select Table Browser under Tools
    - Configure Assembly hg19, Table knownCanonical, Output filename hg19-kc.tsv
    - Use `mv` to move the file from `~/Desktop` to `~/qbXX-answers/weekX`
    - Confirm that `hg19-kc.tsv` has 80,270 lines
    - Use the following command to convert from a .tsv to a .bed file
        ```
        cut -f1-3,5 hg19-kc.tsv > hg19-kc.bed
        ```
    - Confirm that the first few lines of `hg19-kc.bed` looks like
        ```
        #chrom	chromStart	chromEnd	transcript
        chr1	169818771	169863037	ENST00000367771.11_11
        chr1	169764180	169823221	ENST00000359326.9_7
        chr1	27938574	27961696	ENST00000374005.8_7
        ```

    Count how many genes are in each 1 mb interval using `bedtools intersect`

    - You will need to use exactly one 1-letter option e.g.
        ```
        bedtools intersect -x -a fileA.bed -b fileB.bed
        ```
    - Save the output to `hg19-kc-count.bed`
    - Confirm that hg19-kc-count.bed has 3,114 lines

    Use R to plot gene density across all 25 main chromosomes

    - Load `tidyverse`
    - Read in `hg19-kc-count.bed` using `read_tsv()`, specifying the header e.g.
        ```
        header <- c( "chr", "start", "end", "count" )
        df_kc <- read_tsv( ________, col_names=header )
        ```
    - Plot using `ggplot()`, `geom_line()`, and `facet_wrap()`, using the `scales` argument to let the x- and y-axes adjust for each chromosome
        - If you need a place to start, plot just chr1 by first subsetting using `filter( chr == "chr1" )` so that you don't have to use `facet_wrap()`, then proceed to plotting all the chromosomes
    - Save your plot using `ggsave( "exercise1.png" )`

    Submit the following

    - hg19-main.chrom.sizes, hg19-1mb.bed, and hg19-kc-count.bed
        - Do **not** submit the `hg19-kc` files as they are large
    - exercise1.sh -- Bash script with your makewindows and intersect command
    - exercise1.Rmd -- R Notebok with your code
    - exercise1.png -- Plot of gene density across hg19

1. Compare hg19 gene annotations with hg16 [[fn. 1](#footnotes)]

    Prepare hg16 files as in Exercise 1 with the following modifications

    - Download `hg16.chrom.sizes`
    - Do not rename chrM, only exclude secondary contigs e.g.
        ```
        grep -v _ hg16.chrom.sizes > hg16-main.chrom.sizes
        ```
    
    Visualize both hg19 and hg16 gene distributions on the same line plots

    - Create a combined dataframe using `bind_rows( hg19=dfA, hg16=dfB, .id="assembly" )`
    - Place both distributions on the same plot using `aes( color=___ )`

    Calculate how many genes are unique to each assembly

    - How many genes are in hg19?
    - How many genes are in hg19 but not in hg16?
        - Use `intersect` with a 1-letter option to find genes with no overlaps
    - Why are some genes in hg19 but not in hg16?
    - Answer the same three questions but with respect to hg16

    Submit the following

    - exercise2.sh -- Bash script with your commands to calculate unique genes; place output and answers as comments
    - exercise2.Rmd -- R Notebook with your code
    - exercise2.png -- Plot comparing gene density in hg19 and hg16

1. Explore chromatin states between conditions

1. Identify where variation occurs

## Grading

- Plot gene density -- 3 pt
- Compare hg19 vs hg16 -- 2 pt
- Explore chromatin states -- 3 pt
- Identify variation -- 2 pt

## Footnotes

1. Comparing assemblies requires tools like https://genome.ucsc.edu/cgi-bin/hgLiftOver to account for changes in coordinates.  We ignore this important step to simplify this exercise.
