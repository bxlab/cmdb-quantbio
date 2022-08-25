# Day 2 Lunch: Extend and annotate a BED parser

## Overview

Learning objectives

  - Understand and implement defensive coding practices
  - Clearly identify what each part of your code is doing and why
  - Create a BED parser capable of handling bed12 files
  - Utilize your bed parser in a separate script to perform analysis

## Instructions

  1. Prepare day2-lunch/README.md

    a. Create a `/Users/cmdb/qbb2022-answers/day2-lunch` directory
    b. Use TextMate.app to create a README.md file in the new `day2-lunch` directory
    c. Add the following line plus your answer to excercise 2 to `README.md`:
 
      ```
      # QBB2022 - Day 2 - Lunch Excercises Submission
      ```

  2. Make a working copy of the test data `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/extend_bed_parser/extra_data/hg38_gencodev41_chr21.bed` into your `day2-lunch` directory.
  3. Submit `README.md`, your bed parser script, and your analysis script by pushing them to your `QBB2022-answers` repo on [github.com](http://www.github.com). Submit each part as you finish it rather than waiting until the end. **DO NOT** git add the bed file.

    ```
    git add FILENAME
    git commit -m "My commit message"
    git push
    ```

  4. Verify that the file appears on [github.com](https://www.github.com)

## Excercises

  1. Extend your BED parser

    Create a copy of the script we produced during the live-coding lecture, `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/extend_bed_parser/slides_asynchronous_or_livecoding_resources/bed_parser.py`. This will serve as your starting point for the assignment. Given what you learned in the lecture and the information provided in the [BED documentation](https://samtools.github.io/hts-specs/BEDv1.pdf), create a bed parser capable of handling bed3-bed9 and bed12-formatted files. Your script should include the following features:

      - Fully commented
      - All fields of the correct data type for each entry
      - Parser is a single function that can be imported in another script
      - Report the number of incorrect entries, if any, using `sys.stderr` when loading is complete

    For correctly typing the data, you should make sure that comma-separated lines are converted into lists with entries of the correct data type. Make special note of the fact that UCSC adds an extra comma on the blockSizes and blockStarts lists which needs to be stripped. Your parser should not allow bed10 or bed11 entries, as these are forbidden. It should also make sure that the lengths of blockSizes and blockStarts match blockCount.

    For reporting the number of malformed entries, report this only if the number is greater than zero. This will replace printing a line for each malformed entry. Record the number of malformed records in `README.md`.

  2. Analyze a gene BED file

    You will need to create a second script to perform an analysis of a gene BED file to look at exon counts. Your script should do the following:

      - Import your BED parser function
      - Load and store a bed file
      - Find the median number of exons for genes in the BED file you copied under the instructions

    Here's a hint: for each gene the exons are displayed as blocks. In order to import things from a script (as opposed to a module/library), you just need to reference the script's name without the file extension and the script needs to be in the same directory as the file importing it. Assuming that your parser script was named `myBedParser.py`, you would import it as follows:

    ```
    import myBedParser
    myBedParser.bed_parser()
    ```

    or

    ```
    from myBedParser import bed_parser
    bed_parser()
    ```

  Record the median number of exons in `README.md`.

