# QBB2022 - Day 1 - Lunch Exercises

## Overview

Biological Learning Objectives
- Understand how .bed files store annotations
- Understand how .vcf files store variants
- Explore common human genetic variation

Computational Learning Objectives
- Interact with files using the command line
- Explore text files using Unix programs
- Document your work on GitHub

## Instructions

Add each of your answers to `/Users/cmdb/qbb2022-answers/day1-lunch/README.md`.  Please `git push` after each update and do not wait until the end of the session e.g.

```
# Save answer using TextMate
git add README.md
git commit -m "Add answer for exercise 1"
git push
# Confirm at github.com
``` 

Please run the `github_script.sh` if you have not already.  This is required to configure git on your laptop and link your `qbb2022-answers` directory with a similarly named repository at github.com.  Confirm that your home directory contains the following:

```
/Users/cmdb
├── cmdb-quantbio
├── data
│   ├── bed_files
│   ├── vcf_files
├── qbb2022-answers
```

## Exercises

1. Prepare day1-lunch/README.md

    a. Create a `/Users/cmdb/qbb2022-answers/day1-lunch` directory.

    b. Use TextMate.app to create a `README.md` file in the new `day1-lunch` directory.

    c. Add the following lines and fill in your answer to `README.md`:

    ```
    # QBB2022 - Day 1 - Lunch Exercises Submission

    1. I’m excited to learn <your_short_answer>.
    ```

    d. Submit answer via:

    ```
    git add README.md
    git commit -m "Add README.md for exercise 1"
    git push
    ```

    e. Check submission on https://github.com.

2. Calculate the average number of exons per gene

    a. Make a working copy of `genes.chr21.bed` and `exons.chr21.bed` from `~/data/bed_files` in `day1-lunch`.

    b. What is the mean number of exons per gene?  Add your answer and the commands you used to `README.md`.

    c. Describe briefly what you would do to find the median.

3. Tally chromHMM states in E116 lymphoblastoid cells

    a. Make a working copy of `chromHMM.E116_15_coreMarks_hg38lift_stateno.chr21.bed`.

    b. How many regions of the genome are classified for each state?  Add your tally and the commands you used.

    c. Describe briefly how you would determine which state comprises the largest fraction of the genome.  HINT: consider all the columns.

4. Tally populations among 1000 Genomes Project samples

    a. Make a working copy of `integrated_call_samples.panel`.

    b. How many samples are in each of the pop(ulations) of the AFR super_pop(ulation)?  Add your answer and the commands you used.

    c. Describe briefly how you would do this for all five populations.

5. Explore SNP allele frequencies

    a. Make a working copy of `random_snippet.vcf`

    b. Create a `HG00100.vcf`file (remove HG00096-99, HG00101-NA21144).

    c. How many `0|0`, `0|1`, `1|0`, and `1|1` values are present for HG00100?

    d. How many rows contain AF=1?

    e. How many times can AF=1 appear per row?

    f. Describe briefly how you would extract the AFR values.

## Just for fun

A. What is the most common SNP (e.g. A->C)?

B. Use `cut` to extract just the AFR_AF field from column 8.

C. How many unique AFR_AF values are there?

D. How many unique values are there for the NS= field?

