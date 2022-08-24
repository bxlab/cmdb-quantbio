# Assignment 1: Genome Assembly
Assignment Date: Friday, Sept. 9, 2022 <br>
Due Date: Thursday, Sept. 15, 2022 @ 11:59pm <br>

## Lecture

Slides are available here: [GenomeAssembly.pdf](https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/genome_assembly/slides_asynchronous_or_livecoding_resources/GenomeAssembly.pdf)

## Assignment Overview

In this assignment, you are given a set of unassembled reads from a mysterious pathogen that contains a
secret message encoded someplace in the genome. The secret message will be recognizable as a novel insertion of sequence not found in the reference. Your task is to assess the quality of the reads, assemble the genome, identify, and decode the secret message. If all goes well the secret message should decode into a recognizable english text, otherwise doublecheck your coordinates and try again.

For tips on how to run the tools that you need for this assignment, look in the `Resources` section below.

Finally, keep track of all commands that you use in a text or markdown file. You will need to submit this along with any graphs or images.

## Data

Download the reads and reference genome from:
[https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/genome_assembly/extra_data/asm.tgz](https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/genome_assembly/extra_data/asm.tgz). You can click the link, or you can use `wget`.

After downloading the zipped folder, you'll need to extract it with `tar -zxvf <filename.tgz>`

Note we have provided both paired-end and mate-pairs reads (see included README for details).
Make sure to look at all of the reads for the coverage analysis and kmer analysis, as well as in the assembly.

## Assignment

#### Question 1. Coverage Analysis

- Question 1a. How long is the reference genome? [Hint: Try `samtools faidx`]
- Question 1b. How many reads are provided and how long are they? Make sure to measure each file separately [Hint: Try `FastQC`]
- Question 1c. How much coverage do you expect to have? [Hint: A little arthmetic]
- Question 1d. Plot the average quality value across the length of the reads [We want a screenshot from `FastQC`]

#### Question 2. Kmer Analysis

Use `Jellyfish` to count the 21-mers in the reads data. Make sure to use the "-C" flag to count cannonical kmers,
otherwise your analysis will not correctly account for the fact that your reads come from either strand of DNA.

- Question 2a. How many kmers occur exactly 50 times? [Hint: try `jellyfish histo`]
- Question 2b. What are the top 10 most frequently occurring kmers [Hint: try `jellyfish dump` along with `sort` and `head`]
- Question 2c. What is the estimated genome size based on the kmer frequencies? [Hint: upload the jellyfish histogram to [GenomeScope](http://genomescope.org) and report the min "Genome Haploid Length" in the "Results" section]
- Question 2d. How well does the GenomeScope genome size estimate compare to the reference genome? [Hint: In a sentence or two]

#### Question 3. De novo assembly

Assemble the reads using `Spades`. Spades will *not* run on Windows you must use a linux or mac environment.

- Question 3a. How many contigs were produced? [Hint: try `grep -c '>' contigs.fasta`]
- Question 3b. What is the total length of the contigs? [Hint: try `samtools faidx`, plus a short script if necessary]
- Question 3c. What is the size of your largest contig? [Hint: check `samtools faidx` plus `sort -n`]
- Question 3d. What is the contig N50 size? [Hint: Write a short script if necessary]

#### Question 4. Whole Genome Alignment

Use `MUMmer` for whole genome alignment.

- Question 4a. What is the average identify of your assembly compared to the reference? [Hint: try `dnadiff`]
- Question 4b. What is the length of the longest alignment [Hint: try `nucmer` and `show-coords`]
- Question 4c. How many insertions and deletions are in the assembly? [Hint: try `dnadiff`]

#### Question 5. Decoding the insertion

- Question 5a. What is the position of the insertion in your assembly? Provide the corresponding position in the reference. [Hint: try `show-coords`]
- Question 5b. How long is the novel insertion? [Hint: try `show-coords`]
- Question 5c. What is the DNA sequence of the encoded message? [Hint: try `samtools faidx` to extract the insertion]
- Question 5d. What is the secret message? [Hint: Run the provided script `dna-decode.py` to decode the string from 5c.]


## Submission

The solutions to the above questions should be submitted as a markdown or text file on Github, in your `qbb2022-answers` repo. For each question, label each subproblem, and include both the exact commands you used, as well as the actual answer. Submit any requested figures as well. 


## Resources

#### [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Raw read quality assessment

```
$ fastqc /path/to/reads.fq
```

#### [Jellyfish](http://www.genome.umd.edu/jellyfish.html) - Fast Kmer Counting

```
$ jellyfish count -m 21 -C -s 1000000 /path/to/reads*.fq
$ jellyfish histo mer_counts.jf > reads.histo
```

#### [GenomeScope](http://www.genomescope.org/) - Analyze Kmer Profile to determine genome size and other properties

GenomeScope is a web-based tool so there is nothing to install. Hooray! Just make sure to use the `-C` when running jellyfish count so that the reads are correctly processed.

####  [Spades](http://cab.spbu.ru/software/spades/) - Short Read Assembler.

Normally spades would try several values of k and merge the results together, but here we will force it to just use k=31 to save time. The assembly should take a few minutes.

```
$ spades.py --pe1-1 frag180.1.fq --pe1-2 frag180.2.fq --mp1-1 jump2k.1.fq --mp1-2 jump2k.2.fq -o asm -t 4 -k 31
```

#### [MUMmer](http://mummer.sourceforge.net/) - Whole Genome Alignment

```
$ dnadiff /path/to/ref.fa /path/to/qry.fa
$ nucmer /path/to/ref.fa /path/to/qry.fa
$ show-coords out.delta
```

**WARNING: nucmer and related tools do not like it if/when you have spaces or special characters ('@') in the path to the binaries***

#### [SAMTools](http://www.htslib.org/) - Extract part of a genome sequence using 'samtools faidx' (this will extract from contig_id bases 1234 through 5678)

```
$ ./samtools faidx /path/to/genome.fa contig_id:1234-5678
```
