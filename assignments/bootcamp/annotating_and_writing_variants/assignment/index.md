# Day 2 Homework: Annotating and Writing Variants

## Overview

Learning objectives

  - Clearly identify what each part of your code is doing and why
  - Utilize your VCF parser in a separate script to load VCF records
  - Create a mapping dictionary from a file of source-target pairs
  - Write a formatted text output file 

## Instructions

1. Prepare day2-homework/README.md

    a. Create a `/Users/cmdb/qbb2022-answers/day2-homework` directory

    b. Use TextMate.app to create a README.md file in the new `day2-homework` directory
    
    c. Add the following line plus your answer to excercise 2 to `README.md`:

    ```
    # QBB2022 - Day 2 - Homework Excercises Submission
    ```

2. Make a working copy of the test data `/Users/cmdb/data/vcf_files/random_snippet.vcf` into your `day2-homework` directory.

3. Make a working copy of the dbSnp data `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/annotating_and_writing_variants/extra_data/dbSNP_snippet.vcf` into your `day2-homework` directory.

4. Submit `README.md`, your commented VCF parser script, your analysis script, and your annotated VCF file by pushing them to your `QBB2022-answers` repo on [github.com](http://www.github.com). Submit each part as you finish it rather than waiting until the end. **DO NOT** git add the original VCF file. In order to git add your annotated VCF file, you will need to add `--force` before the file name since VCF files are listed in `.gitignore` as files not to commit.

    ```
    git add FILENAME
    git commit -m "My commit message"
    git push
    ```

5. Verify that the file appears on [github.com](https://www.github.com)

## Excercises

1. Comment the VCF parser

    a. Create a copy of the script we produced during the live-coding lecture, `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/annotating_and_writing_variants/slides_asynchronous_or_livecoding_resources/vcfParser.py`

    b. Comment each block of code and explain what it is doing and why

2. Annotate a VCF file with SNP labels

    Create an analysis script that will import your VCF parser, load both the 1KGP genomes (`random_snippet.vcf`) and dbSNP (`dbSNP_snippet.vcf`) VCF files, and then annotate the 1KGP variants with the IDs from the dbSNP file, when appropriate. Your script should include the following features:

    - Fully commented
    - Loads the records from `dbSNP_snippet.vcf` and `random_snippet.vcf`
    - Creates a dictionary of positions and IDs from the dbSNP file
    - Replaces the ID in each record from `random_snippet.vcf` with the correct label, if it exists, from your dbSNP dictionary
    - Finds and reports the number of `random_snippet.vcf` records that do not have a corresponding ID in your dbSNP ID dictionary

    Record the number of unlabeled records in `README.md`.

3. Write the first 100 **labeled** records from `random_snipett.vcf` to a file

    Use your parsed and annotated `random_snippet.vcf` records list for this part, **not the original file** (i.e. don't pull correctly formatted lines from the original file for your output). Your code should have the following features:

    - Writes only a single header line with the column labels (i.e. the last line of the original header).
    - Correctly formats each record (matching the format of the original file)
    - Writes only the first 100 lines that have IDs taken from `dbSNP_snippet.vcf`

    You can choose how your code writes to a file, either by printing to sys.stdout and using `>` to put the output into a file or using `open(fname, 'w')` to write to a file within your script. If you use the first one, you can still have your script report the number of unannotated records using `sys.stderr`.
    