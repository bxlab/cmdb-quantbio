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

3. Submit `README.md`, your commented VCF parser script, your analysis script, and your annotated VCF file by pushing them to your `QBB2022-answers` repo on [github.com](http://www.github.com). Submit each part as you finish it rather than waiting until the end. **DO NOT** git add the original VCF file. In order to git add your annotated VCF file, you will need to add `--force` before the file name since VCF files are listed in `.gitignore` as files not to commit.

    ```
    git add FILENAME
    git commit -m "My commit message"
    git push
    ```

4. Verify that the file appears on [github.com](https://www.github.com)

## Excercises

1. Comment the VCF parser
    Identify what each block of code is doing and why.

2. Annotate a VCF file with SNP labels

    Create a copy of the script we produced during the live-coding lecture, `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/annotating_and_writing_variants/slides_asynchronous_or_livecoding_resources/vcf_parser.py`. Create a script that uses the labels from the `dbSnp_snippet.vcf` file and applies them to the correct SNP records in `random_snippet.vcf`. Your script should include the following features:

    - Fully commented
    - Loads the records from `dbSnp_snippet.vcf` and creates a dictionary of positions and labels
    - Labels all records from `random_snippet.vcf` with correct label, if it exists

3. Report the number of records from `random_snippet.vcf` that don't have a label

    Record the number of unlabeled records in `README.md`.

4. Write the first 100 **labeled** records from `random_snipett.vcf` to a file

    Use your parsed vcf records list for this part, **not the original file** (i.e. don't pull correctly formatted lines from the original file for your output). Write only a single header line with the column labels (i.e. the last line of the original header). Also, make sure that each record is formatted correctly (matching the format of the original file). Include this functionality in your annotation script (excercise 2).
