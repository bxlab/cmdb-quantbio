# QBB2024 - Day 2 - Lunch Exercises

## Overview

Biological Learning Objectives
- Summarize gene metadata retreived from [BioMart](https://www.ensembl.org/info/data/biomart)
- Explore the [GENCODE](https://www.gencodegenes.org) genome annotation

Computational Learning Objectives
- Practice working at the [command line](https://swcarpentry.github.io/shell-novice/04-pipefilter.html)
- Explore text files using [core Unix programs](https://en.wikipedia.org/wiki/List_of_POSIX_commands)

## Instructions

Complete the following exercises using a combination of `cut`, `sort`, `uniq`, and `grep`.
Document each command you use, along with any short answers, in a single `README.md` file stored in `~/qbb2024-answers/day2-lunch`.
Use [Markdown](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) headings, lists, and code blocks to organize your answers into sections and format content for proper display e.g.

    ## Answer 1
    
    - `wc -l hg38-gene-metadata-feature.tsv`
    - There are 61633 lines

Please `git push` after each exercise and **do not wait** until the end of the session.

## Exercises

[BioMart](https://www.ensembl.org/info/data/biomart) provides a way to obtain genome annotation information (referred to as Attributes) such as gene IDs, transcripts, phenotypes, GO terms, structures, orthologues, variants, and more.
This information can be retrieved through the [web-based tool](https://www.ensembl.org/biomart/martview) as well as the Bioconductor [biomaRt](https://www.bioconductor.org/packages/biomaRt) package.

1. Tally the number of each `gene_biotype` in hg38-gene-metadata-**feature**.tsv.  How many `protein_coding` genes are there?  Pick one biotype you would want to learn more about and explain why.

2. Which `ensembl_gene_id` in hg38-gene-metadata-**go**.tsv has the most `go_ids`?  Create a new file that only contains rows corresponding to that gene_id, sorting the rows according to the `name_1006` column.  Describe what you think this gene does based on the GO terms.

    ```
    ENSG000XYZ  GO:0090425  acinar cell differentiation
    ENSG000XYZ  GO:0016323  basolateral plasma membrane
    ENSG000XYZ  GO:0045296  cadherin binding
    ```

[GENCODE](https://www.gencodegenes.org) works to annotate the human and mouse genome using biological evidence such as long-read RNA-seq, Ribo-seq, and other targeted approaches.
This gene set is used by many projects including Genotype-Tissue Expression ([GTEx](https://gtexportal.org)), The Cancer Genome Atlas ([TCGA](https://www.cancer.gov/ccg/research/genome-sequencing/tcga)), and the Human Cell Atlas ([HCA](https://www.humancellatlas.org)).

3. Immunoglobin (Ig) genes are present in over 200 copies throughout the human genome.  How many [IG genes](https://www.gencodegenes.org/pages/biotypes.html) (not pseudogenes) are present on each chromosome?  You can use a dot (`.`) in a [regular expression](https://www.regular-expressions.info/quickstart.html) pattern to match any single character.  How does this compare with the distribution of IG pseudogenes?

4. Why is `grep pseudogene gene.gtf` not an effective way to identify lines where the `gene_type` [key-value pair](https://www.gencodegenes.org/pages/data_format.html) is a pseudogene (hint: look for overlaps_pseudogene)?  What would be a better pattern?  Describe it in words if you are having trouble with the regular expression.

5. Convert the annotation from [.gtf format](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) to [.bed format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).  Specifically, print out just the chromosome, start, stop, and gene_name.  As `cut` splits lines into fields based on the tab character, first use [`sed`](https://www.gnu.org/software/sed/manual/sed.html) to create a new file where spaces are replaced with tabs.

    ```
    sed "s/ /\t/g" gene.gtf > gene-tabs.gtf
    ```

## Just for fun

A. Explore the GENCODE [mouse annotation](https://www.gencodegenes.org/mouse) noting similaries and differences with the human annotation

