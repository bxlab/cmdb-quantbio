# Mini-Project -- Calculate Assembly Metrics

## Overview

Biological Learning Objectives
- Explore basic [assembly metrics](https://wikipedia.org/wiki/N50,_L50,_and_related_statistics)
- Examine impact of [sequencing technologies](https://wikipedia.org/wiki/Sequence_assembly#Technological_advances)

Computational Learning Objectives
- Practice git by developing a project
- Obtain datasets from the internet
- Parse through a .fasta file and calculate metrics

## Instructions

Document your answers in `~/qbXX-answers/miniproject-assembly-metrics`

`git push` after each exercise and **do not wait** until the end of the session e.g.

```
git add README.md
git commit -m "Start documentation for this assignment"
git push
# Confirm at github.com
```

## Exercises

1. Start documenting your work with a README.md file

	- Change your working directory to your answers repository e.g. `cd ~/qb25-answers`
	- Make a new directory `miniproject-assembly-metrics` using `mkdir`
	- Navigate into your new directory using `cd`
	- Create a new README.md file using `touch`
	- Open the file in VS Code using `open README.md -a "visual studio code"
	- Add a title (e.g. `# My Project`) and a short description using your own words
	- Push your README.md to GitHub
	- Confirm that things look correct at github.com

1. List available genome assemblies for *C. remanei*

	- Navigate using a web browser to https://parasite.wormbase.org/species.html
	- Scroll down to *Caenorhabditis remanei*
	- Click on each of the four assemblies and copy the URLs for the `Genomic Sequence (FASTA)`
	- Add each URL to your README as a list item (e.g. `- https://ftp.ebi.ac.uk/genome.fa.gz`)
	- Push your README.md to GitHub and confirm that things look correct

1. Create a Bash script to download the genome assemblies

	- Create a new file named getGenomes.sh
	- Add a new line for each assembly in the format `wget URL`
	- Make the file executable using `chmod +x`
	- Run your script using `./getGenomes.sh`
	- Confirm that you have four files that are 30-40 Mb using `ls -lh`
	- Uncompress the four files using `gunzip *.gz`
	- Update your README.md with the **uncompressed** file sizes
	- Push your script and README.md and to GitHub

1. Create a Python script that reads in a single .fasta file and reports basic metrics

	Fetch `fasta` module which defines `FASTAReader` class

	- Download [fasta.py](https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/resources/code/fasta.py) to your directory using `wget`
	- Open fasta.py in VS Code to skim over the FASTAReader class

	Write a Python script that expects a .fasta file as an argument (e.g. `./assembly-metrics.py sequences.fa`)

	- Create a new file named `assembly-metrics.py`
	- Add a "hash-bang" (`#!`) line to declare this a python3 script
	- Import the sys and fasta module e.g. `import fasta`
	- Use `open()` to open the file specified as the "1-th" command line argument e.g. `[1]`
	- Create a FASTAReader object using `fasta.FASTAReader( ___ )`
	- Iterate through FASTAReader using `for ident, sequence in ____`
		- Count the number of contigs
		- Determine each sequence length using `len()` and sum up the total length
	- Print the "Number of contigs: ", "Total length: ", and "Average length: "

	Analyze the four .fa files using your script

	- Run your script on all four .fa files
	- Update README.md with instructions on using your script and summary of results
	- Push your script and updated README.md

1. Extend your Python script to calculate the N50 statistic

	- Create a list to store contig lengths
	- Append each contig length to your list using `___.append( ___ )`
	- Sort contig lengths using `___.sort( reverse=True )`
	- Iterate through the list using a for loop and at each iteration
		- Sum up the cumulative length thus far
		- Stop when the cumulative length is greater than half the total length
	- Print out the "sequence length of the shortest contig at 50% of the total assembly length"	
	- Run your script on all four .fa files and
	- Update README.md with a summary of results
	- Push your script and updated README.md

## Grading

- README -- 2 pt for having all five details (description, URLs, file sizes, first summary, second summary)
- Bash script -- 1 pt for `getGenomes.sh`
- Python script -- 3 pt for original `assembly-metrics.py`
- N50 statistic -- 2 pt for extending Python script
- git -- 2 pt for making a commit after each of the five exercises

