# Day 3 Homework: Parsing and analyzing GFF files

## Overview

GFF stands for General Feature File, and is a flexible file for describing genomic annotations. In this case, you will be using GFF files that describe pseudogene annotations called by RefSeq.

Learning objectives

  - Understand the format and purpose of a GFF file
  - Parse a GFF file, including converting the "INFO" field into a dictionary of attributes
  - Using two dictionaries to identify common or unique sets of keys in each condition's dictionary
  - Plot results for sets of unique and shared features

## Background

The complete human genome was recently released from the T2T consortium, filling in all the remaining gaps and removing much incorrectly-placed repetitive sequence (chm13v2.0). You will be looking at the impact of this new assembly on the annotation of pseudogenes throughout the genome compared to the latest update from the previously most up to date human genome assembly (GRCh38).

## Instructions

1. Prepare day3-homework/README.md

    a. Create a `/Users/cmdb/qbb2022-answers/day3-homework` directory

    b. Use Sublime.app to create a README.md file in the new `day3-homework` directory

    c. Add the following line plus your answers to excercises 6 and 7 to `README.md`:

    ```
    # QBB2023 - Day 3 - Homework Excercises Submission
    ```

2. Make a working copy of the test data `/Users/cmdb/cmdb_quantbio/assignments/bootcamp/parsing_gff/extra_data/*.gff` into your `day3-homework` directory.

3. To the `QBB2023-answers` repo on [github.com](http://www.github.com), submit the following:
    
    Submit `README.md`
    Your commented GFF parser/plotting script
    Your analysis plot
    Answers to exercises 6 and 7

Submit each part as you finish it rather than waiting until the end. **DO NOT** git add the original GFF file.

```
git add FILENAME
git commit -m "My commit message"
git push
```

4. Verify that the file appears on [github.com](https://www.github.com)

## Exercises

1. Parse the two GFF files

    a. Because each file cshould be identical in terms of it's format, this is an ideal use case for a function to reuse your code

    b. For each line in the file, you will need to remove the newline character at the end of the line *.rstrip("\n")* and split the file at each tab *.split("\t")*. We are only concerned with the pseudogenes, which is specified by the third column (index 2 in the split line, since Python starts counting at zero)

    c. The only field that we really care about once filtering out non-pseudogenes is the INFO field, which is the last field. It contains lots of information crammed into a long string. You will need to break this line up into its own dictionary. Each entry is separated by a semi-colon ";", and the key and value for each entry is separated by an equal sign. You can create this dictionary by first creating a list using *.split(";")*. Then you can use a *for* loop, and splitting each list entry to create the key and variable pair.

    d. Note that this next step is a little tricky so pay attention. There are two additional pieces of information we're interested in. First, the size of the pseudogene (note the start and stop coordinate columns). Second, there are some gene names that appear multiple times in each file. Each gene should have a single entry in your gene dictionary but you want to keep track of how many times it appeared in the original file. The easiest way to handle this is to add a "length" and "count" key to the INFO dictionary that you created by parsing the info field.

    e. Add your newly parsed INFO dictionary into your master gene dictionary. You will be using the gene name as the key (in the info field, there is a key: value pair of "gene": NAME). Before adding each pseudogene entry into your master dictionary of pseudogenes, you need to check if there is already an entry for that given gene. In cases where the gene already has an entry, you should simply add one to that gene's count. Also, let's only keep the longest length seen for a given gene name.
    ```
    fs = open(fname, 'r')
    genes = {}
    for line in fs:
        ...
        key = info['gene']
        if key in genes:
            genes[key]['count'] += 1
            genes[key]['length'] = max(genes[key]['length'], info['length'])
        else:
            genes[key] = info
    fs.close()
    ```
    
2. Determine which genes are shared between the two annotations, and which are unique to each annotation set. This is a good time to use a *for* loop with *dict.keys()* or *dict.items()* along the the *in* conditional test.

3. For each annotation and the set of shared genes, find genes that occurred more than once in the original annotation file (i.e. they appear to have more than one genomic copy). For each one, keep track of the copy number and length of the gene.

4. For each annotation and its set of unique genes, keep track of the copy number and length of each gene.

5. Create a 4x4 set of panels in a single plot. In the first two panels, plot the shared multi-copy genes with copy number on the x-axis and length on the y-axis. In the last two panels, plot the unique genes with copy number on the x-axis and length on the y-axis. Make sure to add a descriptive title to each panel and preferably label the x and y axes. Because you will be using the same plotting approach for each, this is again a great time to create a function to plot each panel.
    ```
    # To work with a multi-panel plot, you can create it as follows:
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    # To plot to a specific panel, we need to index ax
    # with the first axis specifying row and the second specifying column
    ax[0, 0].scatter(Xdata, Ydata, alpha=0.2)
    # Use alpha=0.2 so be can see overlapping points better

    ax[0, 0].set_title("A descriptive title")
    ...
    plt.tight_layout() # This optimizes the margins for prettier plots
    plt.figsave("FIGURENAME.png")
    ```
6. Report how many genes had multiple putative copies and what the maximum number of copies was for each annotation.

7. Look at the plots and try to make an reasonable guess as to what they tell you about the effects of resolving highly repetitive and unplaced sequences in chm13v2.0 vs. GRCh38.

***

## Advanced Exercises

1. Make all of your plots have the same x-axis and y-axis limits so that all data are still included and on the same scale between panels.

2. Repeat the above exercise, but further restrict your analysis to only "Curated Genomic" pseudogenes. How did this affect the results?

3. Using functions, get rid of redundant code in your script such that any repeated action is handled by a function.

4. For the unique set and the shared set of genes for each annotation, tally up how many times each gene description appears. Print out all descriptions tied for most occurrences for each of the four conditions (shared and unique for aech annotation).

