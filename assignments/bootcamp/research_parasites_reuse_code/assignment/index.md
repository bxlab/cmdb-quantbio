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

1. Reproduce plots using https://github.com/bxlab/cmdb-plot-vcfs

    Add the following to your `day4-lunch/README.md` file:

    - Portion of `do_all.sh` output that reports how many bp each feature covers
    - Strategy to confirm that reproduced plots are the same as in the `cache/` directory
    - Three other `gene_type`s present in the GENCODE .gtf that you find interesting and why 

2. Modify workflow

    Within the `cmdb-plot-vcfs/` directory:

    - Improve individual plots e.g. log scale, same y-axis, title
    - Create lncRNA.chr21.bed.vcf.png plot

    Within `day4-lunch/` directory:

    - Commit the improved plots and files you modified
    - Describe possible trends among plots in `README.md`

3. Create documentation

    Add documentation for bxlab/cmdb-plot-vcfs in `day4-lunch/README.md` including:

    - Synopsis -- <50 words
    - Usage -- syntax including input file requirements
    - Dependencies -- software requirements
    - Description -- how it works (bullet points or prose)
    - Output -- example output

## Just for fun

A. Modify workflow to create a single plot e.g. multiple series or subpanels

B. Browse the commit history at github.com and explain what was wrong with the initial bp calculation.  Checkout the original calculation, re-run the workflow, and document the original bp coverage caluclations.

C. Finish https://github.com/git-game/git-game (and also https://github.com/git-game/git-game-v2)
