# QBB2022 -- Day 1 Homework: Genome Arithmetic

## Prepare `day1-homework` directory

a. Create a `/Users/cmdb/qbb2022-answers/day1-homework` directory.

b. Make working copies of all of the scripts included in `/Users/cmdb/cmdb-quantbio/assignments/bootcamp/bedtools_genome_arithmetic/scripts` within your new `day1-homework` directory.

c. Use TextMate.app to create a `README.md` file in the new `day1-homework` directory.

d. Add the following line to `README.md` and save the file. Record your answers to exercises throughout the homework assignment within README.md

```
# QBB2022 - Day 1 - Homework Exercises Submission
```
e. Submit your answers to your `qbb2022-answers` repository as you work on the assignment and complete an exercise; do not wait until all of the exercises are complete. Use informative commit messages. For example:

```
git add README.md exercise1.sh
git commit -m "edited script and answers for day 1 hw exercise 1"
git push
```
f. Verify that your files were successfully uploaded to your online repository on GitHub.  

## Exercise 1

In the interactive lecture, we looked at what's the most common alternate allele for a reference allele of Cytosine using `awk`. You are going to use the script `exercise1.sh` in order to find the most common alternate allele for all of the reference allele bases, not just Cytosine. However, there is an error in the script that needs to be debugged.

  * Look at the script and compare it to the command we used in the interactive lecture. Using the usage statement provided at the beginning of the script, try to run the script using `~/data/vcf_files/random_snippet.vcf` as the input VCF file.
  * What error message do you see and how will you fix it? If needed, google part of the error message to see how you can fix the script.
  * Record the output from the working script.
  * Do the results make sense given what you know about the biology of these bases?

## Exercise 2

What's the most common alternate allele for a Cytosine reference allele for variants occurring in promoter like regions of the genome?

  * Use the `~/data/vcf_files/random_snippet.vcf` file as the variant informations for this problem.
  * To define promoter like regions of the genome, use the provided chromHMM segmentation in `~/data/bed_files/chromHMM.E116_15_coreMarks_hg38lift_stateno.chr21.bed` and consult [this resource](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state).
  * Using this segmentation, do promoters appear to be clearly and objectively defined?
  * Use both `bedtools` and `awk` to answer this question. Record any reasoning, as well as your commands and output in your README file.
  * Finally, do the observed results lead you to any biological observations, hypotheses, or conclusions?

## Exercise 3

It is often of interest for biologists to find the closest gene to some variant or set of variants which are identified in a screen or as a significant GWAS hit. Debug the script that uses `bedtools closest` to find the closest gene for each variant.  

  * Look at the script (`exercise3.sh`) and paraphrase in your README what you think each line is doing and why it is necessary given the documentation for `bedtools closest`.
  * Using the usage statement provided at the beginning of the script (`exercise3.sh`), try to run the script using `~/data/vcf_files/random_snippet.vcf` as the input VCF file.
  * There are two errors in the script. Once you successfully debug the first and run the script again, you will see the second error. For each error, record what error message you see and how you will fix it. If needed, google part of the error message to see how you can fix the script.
  * Submit the corrected script.
  * Finally, use the output to answer the following questions. Record the answers and the commands you used to find the answers in your README.
    * How many variants are returned and how many unique genes are returned? How many variants on average are therefore connected to a gene with `bedtools closest`?

## Exercise 4

You've been given a script `exercise4.sh` which intends to first use `bedtools intersect` in order to find the genes that intersect with H3K27ac in naive B cells and the genes that intersect with H3K9me3 in naive B cells. It saves the output of each of these intersections into separate files. Then the script intends to perform a comparison between the two resulting gene/histone mark intersection files, reporting any genes that uniquely intersect with H3K27ac but never intersect with H3K9me3 within naive B cells.

Look at the script and annotate within your README file what you think it is doing at each step.

In this fourth exercise, the script `exercise4.sh` has both syntax errors and a logic error. Fix the script's syntax errors so that it runs and the logic error so that it performs the desired operation.

Using the usage statement provided at the beginning of the script (`exercise4.sh`), try to run the script, record the following in your README, and submit the fixed script:

  * What were the syntax errors and what was the logic error?
  * What is the output of the script when it follows its intended goal

## Exercise 5

In the `~/data/bed_files` directory, you have bedgraph files for 6 different Histone marks (H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3, and H3K9me3) and 3 different cell types (CD4, CD8, and naive B cells).

Use the appropriate files and `bedtools` subcommand to find the regions that are marked by H3K36me3 in CD8 cells, but not by H3K36me3 in CD4 cells and not by H3K36me3 in naive B cells. If there is any overlap, do not consider the CD8 region to be unique.

Within your README, record

  * What fraction of the original H3K36me3 regions in CD8 cells were completely unique to CD8 cells when compared to CD4 and naive B cells?
  * Your `bedtools` command and any commands that you used to answer these questions.


## Advanced Exercise 1 (optional)

Using the bedgraph files in `~/data/bed_files` again, specifically, within just the CD4 cell files, use the appropriate `bedtools` subcommand to see which regions share some or all of these marks. Save your output with an informative header so that you can further explore the results.

Within your REAMDE, record

  * How many regions had some amount of overlap between the 6 histone marks?
  * How many regions were marked by all 6 histone marks?
  * Which combination of marks occurred most frequently?
  * The commands that you used to answer all of these questions

## Advanced Exercise 2 (optional)

In the interactive lecture, we found the average inter-variant distance for biallelic variants on chromosome 21. Using that `awk` command as reference, build and use the appropriate input file to find the average inter-variant distance for biallelic variants on chromosome 21 which occur within human genes.


## Advanced Exercise 3 (optional)

In the interactive lecture, we saw that a handful of bedtools subcommands didn't have an online manual. Pick any one of these and can you find information or resources that describe how the tool works?
