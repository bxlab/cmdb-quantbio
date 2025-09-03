# Unix and Python Scripts

## Overview

Biological Learning Objectives
- Explore genome annotation files in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [GTF format](https://ensembl.org/info/website/upload/gff.html)

Computational Learning Objectives
- Use [Unix commands](https://wikipedia.org/wiki/List_of_POSIX_commands) to examine text files
- Create Python scripts for text wrangling

## Part I -- Unix

### Instructions

Work on these exercises locally on your laptop in `~/qbXX-answers/unix-python-scripts`

Submit your answers via https://github.com

- Create a new README.md in your repository under a new directory named `unix-python-scripts`
- Commit your answer after each exercise and do not wait until the end of the session
- Place each of your Unix commands along with the output as a comment e.g.

    ```
    tail ce11_genes.bed | head -n 2

    # chrIII	13768540	13771741	NM_067444.8	515	-
    # chrIII	13769876	13769953	NR_003432.1	9	-
    ```

### Exercises

1. Explore ce11_genes.bed using Unix

    Calculate each of the following statistics by constructing a single command using one or more (linked together with a `|`) of the following commands: `cut`, `grep`, `sort`, `uniq`, `wc`

    - How many features (lines)?
    - How many features per chr? e.g. `chrI`, `chrII`
    - How many features per strand? e.g. `+`, `-`

1. Explore GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt using Unix

    - Which three SMTSDs (Tissue Site Detail) have the most samples?
    - How many lines have "RNA"?
    - How many lines do **not** have "RNA"?

1. Explore ~/Data/References/hg38/gencode.v46.basic.annotation.gtf using Unix

    - How many entries are there for each feature type?  Look at column 3 and be sure to skip any lines that begin with `#`
    - How many lncRNA entries are on each chromosome?

## Part II -- Python Scripts

### Instructions

Document your answers in `~/qbXX-answers/unix-python-scripts`

Upload your scripts to https://github.com after each exercise and do not wait until the end of the session

### Exercises

1. Recalculate ce11_genes.bed scores using Python

    Write a script that for each feature (line) recalculates the score (column 5) such that

    - new_score = original_score * feature_size
    - new_score is positive or negative based on the strand (column 6)

    Print out all six columns in BED format

1. Export gene features to BED format using Python

    Write a script that takes gencode.v46.basic.annotation.gtf and

    - Ignores a line if it startswith `#`
    - Only prints output for `gene` features (column 3)
    - Subtracts 1 from the start position to convert to [zero-based coordinate system](https://en.wikipedia.org/wiki/BED_(file_format)#Coordinate_system)
    - Correctly parses gene name from the `attribute` (column 9)
    - Prints chr, start, stop, gene name

1. Transform GTEx data using Python

    Write a script that extracts expression values for the first gene (DDX11L1) which is stored on a single line spread across more than 17,000 columns and transposes the data so that the expression in each sample is stored on a separate line.

    Open the expression data file GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct 
    - Skip the first two lines using `.readline()`
    - Read the next header line and save the fields after splitting
    - Read the next data line and save the fields after splitting
    - Create a dictionary by looping through the fields, using `header[i]` as the key to store `data[i]` as the value

    Open the metadata file GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
    - If the SAMPID in the first column is present in the dictionary, print out the SAMPID, expression, and SMTSD e.g.

    ```
    GTEX-ZZPU-2426-SM-5E44I 0       Artery - Tibial
    GTEX-ZZPU-2626-SM-5E45Y 0.01965 Muscle - Skeletal
    GTEX-ZZPU-2726-SM-5NQ8O 0.02522 Adipose - Subcutaneous
    ```

    What are the first three tissues that have >0 expression?

## Grading

- unix-commands.sh -- 3 pt
- recalculate-score.py -- 2 pt
- sample-by-gene.py -- 3 pt
- gtf2bed.py -- 2 pt

