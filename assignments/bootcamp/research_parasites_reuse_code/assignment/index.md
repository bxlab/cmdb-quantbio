# QBB2022 - Day 4 - Lunch Exercises

## Overview

Biological Learning Objectives
- Explore size of chromosome features (e.g. exons, processed_pseudogene)
- Compare allele counts among these regions

Computational Learning Objectives
- Practice reproducible research
- Orchestrate analysis using Bash scripts
- Customize matplotlib plots

## Instructions

Add each of your answers to `qbb2022-answers/day4-lunch`.  Please `git push` after each update and do not wait until the end of the session.

## Exercises

1. Reproduce plots using [bxlab/cmdb-plot-vcfs](https://github.com/bxlab/cmdb-plot-vcfs)

    Add the following to your `day4-lunch/README.md` file:

    - Portion of `do_all.sh` output (not script) that reports how many bp each feature covers
    - One or more strategies to confirm that reproduced plots are the same as in the `cache/` directory
    - Three other `gene_type`s present in the GENCODE .gtf that you find interesting and why 

2. Modify workflow

    Work inside the `cmdb-plot-vcfs/` directory to:

    - Improve individual plots e.g. log scale, same y-axis, title
    - Create lncRNA.chr21.bed.vcf.png plot

    Change into the `day4-lunch/` directory and:

    - Copy over the improved plots and files you modified
    - Commit only these plots and files
    - Describe possible trends among plots in `README.md`

3. Create documentation

    Add documentation for bxlab/cmdb-plot-vcfs in `day4-lunch/README.md` including:

    - Synopsis -- <50 words
    - Usage -- syntax including input file requirements
    - Dependencies -- software requirements
    - Description -- how it works (bullet points or prose)
    - Output -- example output

    Look to other documentation like [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) for inspiration
    
    ```
    SYNOPSIS
        bxlab/cmdb-plot-vcfs -- ...

    USAGE
        bash do_all.sh <thing1> ...

        <thing1>   ...

    DESCRIPTION
        1. Create .bed files for features of interest
            - Run subset_regions.sh Bash script
            - Use grep to ...
    ```

## Just for fun

A. Modify workflow to create a single plot e.g. multiple series or subpanels

B. Browse the commit history at github.com and explain what was wrong with the initial bp calculation.  Checkout the original calculation, re-run the workflow, and document the original bp coverage calculations.

C. Finish [git-game](https://github.com/git-game/git-game) (and also [git-game-v2](https://github.com/git-game/git-game-v2))
