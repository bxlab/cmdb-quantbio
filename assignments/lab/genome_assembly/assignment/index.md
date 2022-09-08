# Assignment 1: Genome Assembly
Assignment Date: Friday, Sept. 9, 2022 <br>
Due Date: Friday, Sept. 16, 2022 @ 1:00pm ET<br>


## Lecture

**Slides** are available here: `~/cmdb-quantbio/assignments/lab/genome_assembly/slides_asynchronous_or_livecoding_resources/GenomeAssembly.pdf`


## Assignment Overview

In this assignment, you will first implement a short programming exercise and then perform a genomic analysis from short read data. For the genomic analysis you are given a set of unassembled reads from a mysterious pathogen that contains a secret message encoded someplace in the genome. The secret message will be recognizable as a novel insertion of sequence not found in the reference. Your task is to assemble the genome, then identify and decode the secret message. If all goes well, the secret message should decode into recognizable english text, otherwise double-check your coordinates and try again.

For tips on how to run the tools that you need for this assignment, look in the `Resources` section below.

Create a directory called `week1-homework` in your `qbb2022-asnwers` directory. Keep track of all commands that you use in a README markdown file. You will need to submit this along with any graphs or images.


## Data

The reads and reference genome are in a zipped folder here: `~/cmdb-quantbio/assignments/lab/genome_assembly/extra_data/asm.tar.gz`. Copy this file to the `week1-homework` directory you made for this assignment.

After copying the zipped folder, you'll need to extract it with `tar -zxvf <filename.tar.gz>`. You should get 7 files:
1. dna-decode.py
2. frag180.1.fq
3. frag180.2.fq
4. jump2k.1.fq
5. jump2k.2.fq
6. README
7. ref.fa

Note we have provided both paired-end and mate-pairs reads (see included README for details).


## Assignment

#### Question 1. Coverage simulator

- **Question 1.1**. How many 100bp reads are needed to sequence a 1Mbp genome to 5x coverage? How many are needed for 15x coverage? 

- **Question 1.2**. Write a program to simulate sequencing 5x coverage of a 1Mbp genome with 100bp reads and plot the histogram of coverage. Note you do not need to actually output the sequences of the reads, you can just randomly sample positions in the genome and record the coverage. You do not need to consider the strand of each read. The start position of each read should have a uniform random probabilty at each possible starting position (1 through 999,901). You can record the coverage in an array of 1M positions. Overlay the histogram with a Poisson distribution with lambda=5

- **Question 1.3**. Using the histogram from Q1.2, how much of the genome has not been sequenced (has 0x coverage)? How well does this match Poisson expectations?

- **Question 1.4**. Now repeat the analysis with 15x coverage: 1. simulate the appropriate number of reads, 2. make a histogram, 3. overlay a Poisson distribution with lambda=15, 4. compute the number of bases with 0x coverage, and 5. evaluate how well it matches the Poisson expectation.<br>

#### Question 2. De novo assembly

Assemble the reads using `Spades`. <!---Spades will *not* run on Windows you must use a linux or mac environment.-->

- **Question 2.1**. How many contigs were produced? [Hint: try `grep -c '>' contigs.fasta`]

- **Question 2.2**. What is the total length of the contigs? [Hint: try `samtools faidx`, plus a short script if necessary]

- **Question 2.3**. What is the size of your largest contig? [Hint: check `samtools faidx` plus `sort -n`]

- **Question 2.4**. What is the contig N50 size? [Hint: Write a short script if necessary]<br>

#### Question 3. Whole Genome Alignment

Use `MUMmer` for whole genome alignment.

- **Question 3.1**. What is the average identify of your assembly compared to the reference? [Hint: try `dnadiff`]

- **Question 3.2**. What is the length of the longest alignment [Hint: try `nucmer` and `show-coords`]

- **Question 3.3**. How many insertions and deletions are in the assembly? [Hint: try `dnadiff`]<br>

#### Question 4. Decoding the insertion

- **Question 4.1**. What is the position of the insertion in your assembly? Provide the corresponding position in the reference. [Hint: try `show-coords`]

- **Question 4.2**. How long is the novel insertion? [Hint: try `show-coords`]

- **Question 4.3**. What is the DNA sequence of the encoded message? [Hint: try `samtools faidx` to extract the insertion]

- **Question 4.4**. What is the secret message? [Hint: Run the provided script `dna-decode.py` to decode the string from 4.3.]


## Submission

The solutions to the above questions should be submitted to Github, in your `week1-homework` directory. You should be submitting three things:
1. The script you wrote for question 1.2 and 1.4
2. A README markdown file containing the commands you ran and the answers for each question (don't copy your script for question 1 into the README). For each question, label each subproblem, and include both the exact commands you used, as well as the actual answer.
3. The two histograms from questions 1.2 and 1.4.

Make sure you **do not** submit any of the files from `asm.tar.gz`. 


## Resources

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
