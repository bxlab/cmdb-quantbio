# Day 2 Lunch: Extend and annotate a BED parser

## Overview

Learning objectives

  - Understand and implement defensive coding practices
  - Clearly identify what each part of your code is doing and why
  - Create a BED parser capable of handling bed12 files
  - Utilize your bed parser in a separate script to perform analysis

## Instructions

1. Prepare day2-lunch

    a. Create a `/Users/cmdb/qbb2022-answers/day2-lunch` directory

    b. Make a working copy of the test data `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/extend_bed_parser/extra_data/hg38_gencodev41_chr21.bed` in your `day2-lunch` directory.

    c. Use TextMate.app to create a README.md file in the new `day2-lunch` directory

    d. Add the following line to your `README.md`:

    ```
    # QBB2022 - Day 2 - Lunch Excercises Submission
    ```

    e. Later, also add your answers from exercise 2 to your `README.md`.

2. Submit the following to your `qbb2022-answers` repo on [github.com](http://www.github.com). Submit the `README.md` everytime you make an addition and each script as you finish it rather than waiting until the end.

    - Extended bed parser script from exercise 1

    - `README.md` with answers to excercises

    - Analysis script from excercise 2

    - **DO NOT** git add the bed file.

    ```
    git add FILENAME
    git commit -m "My commit message"
    git push
    ```

3. Verify that your files appears in your repo on [github.com](https://www.github.com)

## Excercises

1. Extend your BED parser

    Create a copy of the script we produced during the live-coding lecture, `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/extend_bed_parser/slides_asynchronous_or_livecoding_resources/bed_parser.py`. This will serve as your starting point for the assignment. In the interactive lecture, if we could not correctly convert data within a line to the expected data types, then the line was malformed. Now you are going to extend this parser to check for several other conditions in addition to successful data type conversion.

    - Check that the number of fields is appropriate. Your script should allow for files with bed3, bed4, bed5, bed6, bed7, bed8, bed9, and bed12 formatting, but not lines that follow bed10 and bed11 formatting as these entries are prohibited

    - For bed9 and bed12 formatted files, verify that comma-separated entries are converted into lists whose entries are of the correct data type

      - For bed9, verify that there are 3 integers for the itemRGB entries.

      - Make special note of the fact that UCSC adds an extra comma on the blockSizes and blockStarts lists which needs to be stripped

      - For bed12, make sure that the lengths of blockSizes and blockStarts match blockCount

    Make sure your script adheres to the following:

    - Is fully commented

    - The parser remains a single function that can be imported into another script

    - The parser still returns the `bed` list

    - Rather than printing each incorrectly-formatted line, if the number of incorrect or malformed entries is greater than 0, this number is reported using `sys.stderr` after the file has been completely loaded.

    You can find additional format information provided in the [BED documentation](https://samtools.github.io/hts-specs/BEDv1.pdf).

2. Analyze a gene BED file

    Find the median number of exons for genes in the BED file you made a working copy of in the Instructions section. Recall that for each gene, the exons are displayed as blocks, and so the number of exons would be the number of blocks.

    You will need to create a second script to perform this analysis. Your script should do the following:

    - Import your BED parser function

    - Load a BED file

    - Find the median number of exons for genes in the BED file you copied under the instructions

    In order to import things from a script (as opposed to a module/library), you just need to reference the script's name without the file extension and the script needs to be in the same directory as the file importing it. Assuming that your parser script was named `myBedParser.py`, you would import it as follows:

    ```
    import myBedParser
    myBedParser.bed_parser()
    ```

    or

    ```
    from myBedParser import bed_parser
    bed_parser()
    ```

    Record the following in `README.md`.
      - the number of malformed entries within the test data bed file.
      - the median number of exons for genes within the test data bed file.
