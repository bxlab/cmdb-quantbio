# Mini-Project -- Measure Codon Usage

## Overview

Biological Learning Objectives
- Translate [codons](https://en.wikipedia.org/wiki/Genetic_code) into amino acids and tally up each proportion
- Compare [amino acid](https://en.wikipedia.org/wiki/Amino_acid) composition between different types of proteins

Computational Learning Objectives
- Practice working with dictionaries
- More practice developing a project with git

## Instructions

Document your answers in `~/qbXX-answers/miniproject-codon-usage`

`git push` after each exercise and *do not wait* until the end of the session

## Exercises

1. Prepare your new project directory

    - Make a new directory `miniproject-codon-usage`
    - Create a README.md file with a title and short description
    - Use `wget` to download [fasta.py](https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/resources/code/fasta.py), [codons.py](https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/assignments/bootcamp/miniproject-codon-usage/extra-data/codons.py), [cytoplasm.fa](https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/assignments/bootcamp/miniproject-codon-usage/extra-data/cytoplasm.fa), and [membrane.fa](https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/assignments/bootcamp/miniproject-codon-usage/extra-data/membrane.fa)
    - Push **only** your README.md and .py scripts to GitHub, not the two .fasta files

1. Parse CDS sequence into three nucleotide codons

    Write a Python script that expects a .fasta file as an argument (e.g. `./codon-usage.py sequences.fa`)

    - Create a `FASTAReader` based on the `1`-th command line argument
    - Loop through each `ident, sequence` in the .fasta file
    - Now make another loop so that for each `sequence` you
        - Print the `ident`
        - Store the current three nucleotide codon in a variable named `codon`
        - Print the `codon`
        - Advance to the next codon
    - NOTE: You can advance through a sequence either using a `for` loop and position counter or a `while` loop where you repeatedly remove a codon from the sequence until it is empty

    Test your Python scipt using a small subset of sequences

    - Create a test dataset by running `head` on `cytoplasm.fa` and saving the output to `subset.fa`
    - Use this smaller subset as you develop your Python script (e.g. `./codon-usage.py subset.fa`)

    Once you've completed this exercise

    - Comment out the lines that print `ident` and `codon`
    - Add a short paragraph to your README.md describing the output of your script using the subset of sequences, commenting on why you think the output is correct
    - Push your script, `subset.fa`, and README.md to GitHub

1. Extend your code to translate codons into amino acids and count abundance

    Add two dictionaries, one to translate codon to amino acid, one to count amino acid abundance

    - Open codons.py in VS Code to skim over the contents
    - Import the `codons` module which defines a `codons.forward` dictionary
        - Examine this by temporarily adding `print( codons.forward )`
    - Create an empty dictionary named `aas` to count amino acid abundance

    Update your loop to count amino acid abundance

    - Use the `codons` module to translate each codon (key) into amino acid (value)
    - Update the `aas` dictionary to track the occurance (value) of each amino acid (key)
    - Be sure to check if the amino acid exists (e.g.`if key in dictionary`) and if not initialize (e.g. `value = 1`)

    Test your Python scipt using different input sequences

    - Report final counts using a simple `print( aas )`
    - Add a short paragraph to your README.md describing the output of your script using `subset.fa`, commenting on why you think the output is correct
    - Add a short paragraph to your README.md describing the difference in amino acid abundance between `cytoplasm.fa` and `membrane.fa`, commenting on why you think this difference occurs
    - Push your script and README.md to GitHub

1. Improve your code to more easily compare two sets of sequences

    - Accept two .fasta files as arguments (e.g. `./codon-usage.py seqs1.fa seqs2.fa`)
    - Populate two dictionaries that count amino acid abundance (e.g. `aas1` and `aas2`)
        - There are several ways to process two files, from a basic "duplicate your code with adjustments" to more advanced functions and a list of dictionaries
    - Create a tabular report with three columns (e.g. `f"{aa}\t{count1}\t{count2}"`) with rows sorted alphabetically by amino acid.
        - You can created a sorted list of amino acids by starting with `codons.reverse.keys()`, converting to a list with `list()`, and sorting with `sorted()`
    - Save the output comparing `cytoplasm.fa` and `membrane.fa` to a file e.g. `> cyto-v-mem.tsv`
    - Push your script and output to GitHub

1. Examine codon bias

    - Create a new file `codon-bias.py` that expects two arugments, a single .fasta file and an amino acid (single letter, uppercase)
    - Prepare the .fasta file using `FASTAReader`
    - Count codons (not amino acids) using a dictionary named `counts`
    - Report the counts of every codon for the amino acid specified at the command line
        - Review codons.py in VS Code to examine `codons.reverse`
        - Use `codons.reverse` to retrieve the codons (value) for a given amino acid (key)
            - Note that the value is a **list** of codons so you'll need to iterate through each item in the list e.g. value for key `W` is length 1, value for key `L` is length 6
        - Print two columns, with the first being the codon and the second being the count
    - Save the output for `cytoplasm.fa` examining codons for `L` to a file named `bias-L.tsv`
    - Similarly save the output for `A` to `bias-A.tsv`
    - Push your script and output to GitHub

## Grading

- codon-usage.py -- 2 pt for parsing codons, 2 pt for counting abundance, 1 pt for comparing two sequences
- codon-bias.py -- 2 pt for script, 1 pt for output
- git -- 2 pt for making a commit after each of the five exercises
