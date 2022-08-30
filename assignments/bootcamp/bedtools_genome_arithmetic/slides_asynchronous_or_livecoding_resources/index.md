# QBB2022 - Day 1 - Afternoon Interactive Live Coding


## Overview

Computational Learning Objectives

- Practice and expand knowledge of bash command line tools (e.g., cut, grep, uniq, awk) for processing text files
- Introduce programs that can be installed and run from command line (e.g., bedtools) for processing and comparing files
- Further overview and discuss the use of bed and vcf files for storing biological data
- Practice reading manuals and documentation for tools/programs, especially awk and bedtools.
- Practice using Google, StackOverflow, Biostars, StackExchange, GeeksForGeeks for writing and debugging code

Biological and Bioinformatic Analysis Learning Objectives

- Experience the importance of keeping a lab notebook and recording what you did and why
- Understand ways we relate biological data to each other by localizing data to regions in the genome
- Explore annotations of the genome and what they communicate about different genomic regions

Major Resources
- https://www.gnu.org/software/gawk/manual/gawk.html
- https://bedtools.readthedocs.io/en/latest/index.html

Let's discuss syntax errors vs logic errors

## Interactive Exercise 1 -- executable script

1. Within Terminal, let's make a new directory `day1-afternoon` and in there let's write a script (`ie1.sh`).

  ```
  RANDOM=675
  echo $RANDOM "This script ran successfully"
  ```

  What do you expect the output of the script to be when it is run?

2. Within Terminal, we'll try to run the script  (`ie1.sh`) using the following command: `./ie1.sh` or `bash ie1.sh`

  What is displayed in the Terminal and does it match your expectations?<br />
    - If it doesn't match your expectations, is it an error or something else? (If the script was run with `bash ie1.sh`, then you won't see an error, and you may expect to see "675 This script ran successfully" and then be surprised when they see "28128" instead)

3. But `./ie1.sh` didn't work. And what's the error message mean? Let's try to ascertain the problem and fix it so that `ie1.sh` is successfully executed.

  The error message that we see is `-bash: ./ie1.sh: Permission denied`.

  **General Debugging/googling tip:** generalize the error message, removing specifics like your script name or the line where the error occurred (we'll see this in another example later). Use that specific information for yourself.

  **General Debugging/googling tip 2:** Try to paraphrase what you were trying to do. Like in this case we were trying to run a file on command line and permission was denied.

  From googling, we found out that the problem was that the script wasn't executable, so we made it executable by using the command `chmod +x ie1.sh`. Then `./ie1.sh` works.

  **General running a script tip:** use `bash` or `python` in front of a script name if you're not sure if a script will be executable and you know what language the script is run in. Make it executable with `chmod` if you would like to use the `./scriptname` structure.

## Number of SNPs in the 1000 genomes file that intersect each gene and how many unique genes are represented

We're going to write another script. In this script. Let's look at defining variables and calling variables in bash scripts and how that can go wrong.

```
#!/bin/bash

genefile = /Users/cmdb/data/bed_files/genes.bed
echo $genefile
```

And we'll google the error message using a paraphrase

> bash script defining a variable and error command not found

We see that when defining a variable in bash, we can't have spaces around the equal sign.

Let's build out the script some more. We said that we wanted to look at the number of SNPs in the 1000 genomes `random_snippet.vcf` file that are located within a gene.

First, broader overview, we have several annotations for the genome as a whole. We have a gencode annotation file which tells us where genes, and even exons are. It tells us where pseudogenes, long non-coding RNA, etc. We're going to focus on the genes.bed file which is just protein coding genes on chr21.

For a few cell types, we have locations were histone marks were observed. These are located in `~/data/bed_files` and are the bedgraph files

We also have annotations for non-genic functional elements. These are from chromHMM, a model that uses information about histone marks, accessibility, etc. and reports a segmentation of the genome, i.e, an annotation or summary of the expected role a region plays.  

We're going to use bedtools. Let's look at the manual and see what we see and which subcommand we might want to use.

Note there are several subcommands without any sort of info on how to run.

We want `bedtools intersect`

Let's look at the options for `bedtools intersect`. We want to pass a `-a` and a `-b` file, and note that if we just do that, we'll only get info on one file (if the `-a` file is the variant file, we'll only get info on variants, nothing about the genes. We can either switch them, and then we'll get info on the genes only. Or we can add the `-wb` argument, and it'll append the info about the genes on the end.)

```
$vcffile=/Users/cmdb/data/vcf_files/random_snippet.vcf
bedtools intersect -a $vcffile -b $genefile -wb > intersect.out
```
or

```
bedtools intersect -a $genefile -b $vcffile > intersect.out
```

take the `wc -l` of the output to see number of variants.

Which column do we want to look at if we want to see the number of unique genes?

for the original way ... we want to avoid this because I haven't introduced `awk` yet

```
grep -v "#" ../random_snippet.vcf | awk '{print NF; exit}' #2557 this line is telling us the number of fields rather than printing actual output
cut -f 2561 intersect.out | head
```

so then to find number of unique genes...

```
cut -f 2561 intersect.out | sort | uniq | wc -l
```
166 unique genes

if we re-arranged the bedtools command earlier, we would do:

```
cut -f 4 intersect.out | sort | uniq | wc -l
```

What would we do if we wanted to have a list of the unique genes, and a record of the number of times the gene had a variant intersect it?

We'd use `uniq -c` and drop the `wc -l`

## What's the most common alternate allele for a reference allele of C?

In this file, how many of the alternate alleles, had a reference allele that was C?

```
grep -v "#" random_snippet.vcf | cut -f4 | sort | uniq -c
```

>2027 A
2981 C
2931 G
2061 T

now focusing on the cases where the reference allele was C... let's get info about other columns

We're going to use `awk` for this.

What is `awk` and why do we want to use it?

`awk` is fantastic when you need a conditional to subset a file that is setup in a predictable way and we'd like to go through it line by line to extract info. In our case, like making sure that the 4th column in the vcf file, the reference allele, is a specific nucleotide

let's look at the reference sheet a little, point out the basic single quotes, curly braces expression structure. [reference sheet is here](https://bxlab.github.io/cmdb-quantbio/resources/references/unix.html)

We can print all the lines

```
awk '{print}' ~/data/vcf_files/random_snippet.vcf
```

We can print all the lines and ignore the header

```
awk '/^#/{next} {print}' ~/data/vcf_files/random_snippet.vcf
```

This is printing just the fourth column. It's equivalent to using `cut -f 4` on the file, cutting the 4th column from the file.

```
awk '/^#/{next} {print $4}' ~/data/vcf_files/random_snippet.vcf
```

(recall that we can also use the following to skip the header)


```
grep -v "#" ~/data/vcf_files/random_snippet.vcf | awk '{print $4}'
```

Then we can focus on the task at hand which is looking at column 4 (the reference allele) and printing column 5 (the alternate allele) if column 4 has a value we want -- C or Cytosine.

```
awk '/^#/{next} {if ($4 == "C") {print $5}}' ~/data/vcf_files/random_snippet.vcf
```

Finally, we'll summarize or tally these alternate alleles

```
awk '/^#/{next} {if ($4 == "C") {print $5}}' ~/data/vcf_files/random_snippet.vcf | sort | uniq -c
```

> 484 A
 384 G
2113 T

Let's double check/do a sanity check/make sure that this counts every time the reference allele is C.


```
echo $((484 + 384 + 2113))
```

>2981

Let's also put this into a script and see how we can pass command line arguments to a bash script, like you did with python scripts in the prepwork. We're going to pass the vcf file in the script. Command line arguments to a bash script are preceded with a dollar sign and the number used is the index/order of the input argument. Like the first argument is the file name. 

```
#USAGE: bash scriptname.sh input_vcf_file
awk '/^#/{next} {if ($4 == "C") {print $5}}' $1 | sort | uniq -c
```

And let's move the nucleotide of interest to also be an input to the file too in case we'd like to change it easily instead of hardcoding it. It's the second argument after the file name which is the first. so we will use a $2

```
#USAGE: bash scriptname.sh input_vcf_file nucleotide_of_interest

nucoi=$2
awk '/^#/{next} {if ($4 == $nucoi) {print $5}}' $1 | sort | uniq -c
```

**Lecture went through here**

**everything below was not covered in lecture**

Uh oh, we see an error, how can we fix it?

```
nucoi=$2
awk -v nucoi=$nucoi '/^#/{next} {if ($4 == nucoi) {print $5}}' $1 | sort | uniq -c
```

## Number of SNPs that intersect each exon of each gene

You've already seen an exons.bed file in your lunch assignment file, but it only said where exons were, not what genes they were with. Let's remake that file, and then we'll remake that file, but also annotate it with the corresponding gene name.

Let's use `awk` to extract just exons from the gencode annotation file

First, just print the 3 columns.

```
grep "protein_coding" ~/data/gtf_files/gencode.v41.chr21.gtf | awk '{if ($3 == "exon"){print $1,$4,$5}}'
```

Now is the time to make it a bed file. First let's handle the coordinates

```
grep "protein_coding" ~/data/gtf_files/gencode.v41.chr21.gtf | awk '{if ($3 == "exon"){print $1,$4-1,$5}}'
```

Use `cat -t` to see that it's not tab-delimited

```
grep "protein_coding" ~/data/gtf_files/gencode.v41.chr21.gtf | awk '{if ($3 == "exon"){print $1,$4-1,$5}}' | head | cat -t
```

Then we're going to handle the delimiter

```
grep "protein_coding" ~/data/gtf_files/gencode.v41.chr21.gtf | awk 'BEGIN {OFS = "\t"}{if ($3 == "exon"){print $1,$4-1,$5}}' > exons.chr21.bed
```

note that in making this bed file, we did two things, we made sure the output field separator or the delimiter was a tab, and we also made sure to remove one from the start position because bed files are 0 based index and gtf files aren't.

Working toward this: -- looking as we go to make sure we get the genename

Let's take a step back and look at the full annotation file, rather than just the chr21 one. We want exons from protein coding genes, on chr21, and we want to print the same coordinate info we have been printing, but also to get the info on the gene name... so let's look at how we're going to do that...

```
awk 'BEGIN{OFS="\t"}{if ($1=="chr21" && $3=="exon" && $14=="\"protein_coding\";") {gsub(/[";]/,"");print $1,$4-1,$5,$16}}' gencode.v41.annotation.gtf > ~/data/bed_files/exons.chr21.genename.bed
```

Note, if we wanted to have at least only one of those conditions to be true, a mathematical or (`||`) can be used instead of an and (`&&`)

## next to last interactive exercise, using bedtools multiinter and the bedgraph files, if there's time. There may be

Can we find regions that all have one specific histone mark across the 3 provided cell types?

```
bedtools multiinter -i
```


## Interactive Exercise Last if we have time, which we won't

What's the average distance between SNPs in the full VCF file -- let's use `awk`!

How can we do this -- google it?

> distance between SNPs in a vcf file using awk

https://www.biostars.org/p/355034/

Let's understand what it's doing and try to do something similar


```
grep -v "#" ../ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf | awk '{ awkArray[counter++] = $2; } END { for (n=0; n<counter-1;n++) total+=awkArray[n+1] - awkArray[n]} END { print total/NR }'
```

Websites I used

* https://stackoverflow.com/questions/3122442/how-do-i-calculate-the-mean-of-a-column
* https://linuxhint.com/array_awk_command/
* https://www.unix.com/shell-programming-and-scripting/267636-awk-runs-produces-output-but-error.html
* https://www3.physnet.uni-hamburg.de/physnet/Tru64-Unix/HTML/APS32DTE/WKXXXXXX.HTM


```
grep -v "#" ../ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf | awk '{ awkArray[counter++] = $2; } END { for (n=0; n<counter-1;n++) total+=awkArray[n+1] - awkArray[n]} END { print total/NR }'
```

> 42.6605


Double checking our answer with Python

```
>>> import numpy
>>> filename = "../ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf"
>>> locs = []
>>> distance_between = 0
>>> for line in open(filename):
...     if not line.startswith("#"):
...             fields = line.strip('\r\n').split('\t')
...             locs.append(int(fields[1]))
                if len(locs) > 1:
                  distance_between = distance_between + (locs[-1] - locs[-2])

...
>>> distance_between/(len(locs)-1)
42.66059217815314s
```
