---
title: Quantitative Biology Lab Week 13
layout: default
---

# Assignment 10: Metagenomics
Assignment Date: Friday, Dec. 4, 2020 <br>
Due Date: Friday, Dec. 18, 2020 @ 10am ET <br>

## Lecture

[Lecture slides](https://www.dropbox.com/s/q62sjruupybetll/20201204_qblab_metagenomics.pptx?dl=0)

## Homework: Analysis of metagenomic data

Goal: Learn to work with metagenomic sequence data, taxonomically profile communities, and bin assembled fragments into draft genomes. You will be working with a simplified version of a real infant gut microbiome study, which tracked the gut microbiota of a newborn over the course of several days after birth. Please have a glance at the paper [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530670/).

### Step 0: Download and inspect the data

Download the following [tar archive](https://bx.bio.jhu.edu/data/cmdb-lab/week13_data.tar). Hint: Use `tar -xvf` to decompress it.

Inspect the files within the folder. You should see a `READS` folder, containing paired-end read data from 8 samples. For simplicity, you can consider these 8 samples to be gut samples collected from the infant from day 0 to day 8.

You will also see an `assembly.fasta` file, containing the joint assembly of all the samples. Note that entire genomes of the gut bacteria cannot be assembled. Instead, these assembled sequences are the broken up and scrambled pieces of all the various bacterial genomes.

Finally, there is also a `KRAKEN` folder, containing the approximate classifications of both the reads and the assembled scaffolds. KRAKEN is a very fast and powerful tool that uses k-mer profiles instead of direct alignment to match sequences to the database. Normally you would run KRAKEN on the reads and assembly yourself, but since this required 60GB+ of RAM, they were preprocessed for you. These classifications will be useful to us. Use `less -S` to look over the files and make sure you understand what they mean.

### Step 1: Investigate the taxonomic profile of the reads

We would like to view a detailed taxonomic breakdown of all 8 samples to see how the gut microbiota changes across the days following birth. The KRAKEN classifications of the reads should give us this information, but eyeballing the files is not a good idea, so we need a better way to visualize the data.

Think about the best way to do so. You may use any existing tool (such as [Krona-Tools](https://github.com/marbl/Krona/wiki)), or get creative with your own code. You may also need to parse the KRAKEN output to be more usable. (Hint: The `ktImportText` function in the KronaTools package is one way to do this.)

**Submit your best visualizations.**

**Question 1: Briefly comment on the trends you see in the gut microbiota throughout the first week.**

### Step 2: Deconvolute the assembled scaffolds into individual genomes (binning)

Looking at pie charts is a good way to get a feel for what the microbial communities look like, but we want to know more about individual bacteria! In order to track the changes in abundances of individual prokaryotes throughout time, we need to first extract their genomes from the assembly. As mentioned above, their genomes will be highly fragmented, so the challenge is to "bin" these fragments (contigs or scaffolds) into groups. Ideally, each bin will be a single genome. Our grouping needs to be inclusive enough to get all the parts of each genome, but also specific enough to not include the wrong pieces. Sounds hard, right? It is! Luckily, there is software that can help us...

**Question 2: What metrics in the contigs can we use to group them together?**

One software tool that does this is metaBAT. Go ahead and `conda install metabat2`, and look over the usage examples and guidelines [here](https://bitbucket.org/berkeleylab/metabat). Note that you need to align the reads from multiple samples to the assembly to get the coverage of each contig in all the samples. We recommend that you use BWA to align reads to the assembly. Note that this should result in 8 separate alignment files that you will provide to metaBAT.

**Question 3:<br/>(A) How many bins did you get?<br/>(B) Roughly what percentage of the assembly do they represent?<br/>(C) Do you think the sizes of each bin look about right, based on what you know about the size of prokaryotic genomes?<br/>(D) How would you estimate how complete and how contaminated each bin is?**

### Step 3: Estimate the taxonomy of your putative genomes

Rather than using your binning results, please download the desired output of binning [here](https://bx.bio.jhu.edu/data/cmdb-lab/bins.tar). Use this tar archive for Step 3.

Now that you have individual genomes (we hope so, anyway...) we would like to know what they are. Luckily, you already have KRAKEN classifications of all the scaffolds (`assembly.kraken`), so you could use that as a reference. Cross reference the scaffolds in each bin and their respective KRAKEN taxonomies, and come up with your best prediction for what each bin represents.

For example ...
```
$ head bin.1.fa

>NODE_14_length_235766_cov_39.967778
TGAATCACTCTATCTGCTTCTGTTTTTGCTGCTTCAAGTTCATCATGAATTTTAGTCATT
TCATTATTGTATGCAGTCACGCTCGCTGGTTTCTTACCTTCGGTAACTCCCACATGATTT
AATCCTGGACGTTTTTGTTCTAACACTGACTTATCTGCTTTTAATGCATTTTTTGCTTCA
GTTAATTGAGTTAATGCTTCTTCTACTTTCGCTTTCTCATTTCTGATTTCTTCAGCTGTC
GCATCTCCATTGTTAATCACGCCTTGGGCTTCTGCACTAATTCGCTCTGCTTCTACTTTT
TTCGCTCTATAATTATCTGATGTTACTTGTGTCATTCCTGGCGTTGGATCTTGTTCTGCC
GTTGCTTCATCAAGTTGACGTTTTGCTTCAACTAAAGCACTATTATCTGCTTTATTTTGT
AGTAATGAAATTGCATGTTGAATTTTCAATGAGACATCATTAATATTATTTGTCACATTT
GCAACATCAGTATTTGTTGCTCTATCATTTGATAACACTTCATTAGCTTTTCTTCGAACT

$ grep NODE_14_length_235766_cov_39.967778 week13_data/KRAKEN/assembly.kraken

NODE_14_length_235766_cov_39.967778	root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus haemolyticus;Staphylococcus haemolyticus JCSC1435
```

But do this programatically and do something like tally up the various inferred species, family, order, etc.

**Question 4: <br/>(A) What are your predictions for each bin?<br/>(B) This approach to classification is fast, but not very quantitative. Propose one method to more robustly infer the taxonomy of a metagenomic bin.**

### Step 4: Make a heatmap of the individual bin abundances over time

You have your genomes - now we want to see how their abundances shifted in the baby's gut during these 8 days. [Here](https://bx.bio.jhu.edu/data/cmdb-lab/abundance_table.tab) are the abundances of each bin in each sample. This table contains the average bin coverage in each sample.

Use your favorite software package to make a heatmap. All the samples have the same total read count, so you do not need to standardize your abundance values. Label the bins with their corresponding species/genus name (which you determined in Step 3). If you are unsure how to proceed, consider using the seaborn `clustermap` function.

**Submit the heatmap!**

**Question 5: Compare the results of the individual genome abundance analysis to the conclusions you derived from looking at the read taxonomy distributions (from Step 1). Do they agree with each other? What is different?**

## Submit

* All Python code
* Command line record
* Visualizations from Step 1
* Heatmap from Step 4
* Thoughtful answers to Questions 1-5
