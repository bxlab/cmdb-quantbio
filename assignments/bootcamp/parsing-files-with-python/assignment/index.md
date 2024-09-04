# QBB2024 - Day 2 - Homework Exercises

## Overview

Learning Objectives
- Recognize the connection between [Unix commands](https://en.wikipedia.org/wiki/List_of_POSIX_commands) and executable files
- Use [lists](https://docs.python.org/3/tutorial/introduction.html#lists) to store and retrieve data
- Implement logic using [for loops and if statements](https://docs.python.org/3/tutorial/controlflow.html)
- [Parse text files](https://docs.python.org/3/tutorial/inputoutput.html#reading-and-writing-files), manipulate [strings](https://docs.python.org/3/library/stdtypes.html#string-methods), and extract information

## Instructions

Document each answer in a separate .py file stored in `~/qbb2024-answers/day2-homework`.

Please `git push` after each exercise and **do not wait** until the end of the session.

## Exercises

1. `grep.py` -- Warm up by implementing a basic `grep` program to practice indexing the `sys.argv` list, processing a text file line-by-line, and removing the newline character.

    ```
    $ grep.py FIS1 gencode.v46.basic.annotation.gtf
    chr7 HAVANA gene       101239458 101252316 . - . gene_id "ENSG0000021 ...
    chr7 HAVANA transcript 101239472 101245081 . - . gene_id "ENSG0000021 ...
    chr7 HAVANA exon       101244960 101245081 . - . gene_id "ENSG0000021 ...
    ```

2. `gtf2bed.py` -- Create a program that converts genome annotation information from [.gtf format](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) to [.bed format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).  Specifically, print out just the chromosome, start, stop, and gene_name, stripping off both the beginning `gene_name "` and ending `"`.

    ```
    $ gtf2bed.py genes.gtf
    chr1	11869	14409	DDX11L2
    chr1	12010	13670	DDX11L1
    chr1	14696	24886	WASH7P
    ```

3. `tally-fixed.py` -- Starting with `tally.py`, identify and fix the three bugs in this code.  The output of this program should match the output of `cut -f 1 | uniq -c` e.g.

    ```
    $ grep -v "#" gencode.v46.basic.annotation.gtf | cut -f 1 | uniq -c
    197188 chr1
    152253 chr2
    125674 chr3
    ```

## Just for fun

A. Modify your `grep.py` program to accept an optional `-v` flag that inverts the match

- e.g. `grep transcript_id` or `grep -v transcript_id`

B. Implement [tail](https://en.wikipedia.org/wiki/Tail_(Unix))

- Do this **without** loading the entire file into memory using `f.read()`

C. Create a program that for each gene in prints out one line containing the gene_name followed by each transcript_name

```
FIS1 FIS1-201 FIS1-207 FIS1-203
```

D. `cut.py` -- Implement a basic `cut` program where the 1st command line argument specifies the fields to be output, separated by commas .  Use `.split()` to separate this argument into a list, providing you a list that you can loop over to extract the specified fields from each line in the file (i.e. nested `for` loops).  If you store the specified fields in a list, you can combine the entire list into a string using `"\t".join(my_list)`.

    Test this on `hg38-gene-metadata-feature.tsv` to avoid parsing header lines that begin with `#`.

    ```
    $ cut.py 0,2,6 hg38-gene-metadata-feature.tsv
    ensembl_gene_id  chromosome_name  gene_biotype
    ENSG00000228037  1                lncRNA
    ENSG00000142611  1                protein_coding
    ```

