---
title: Quantitative Biology Lab -- Metagenomics
layout: default
---

# Assignment 12: Metagenomics
Assignment Date: Friday, Dec. 9, 2022 <br>
Due Date: Friday, Dec. 16, 2022 @ 1pm ET <br>

## Homework: Analysis of metagenomic data

Goal: Learn to work with metagenomic sequence data, taxonomically profiled communities, and fragments assembled into bins and draft genomes. You will be working with a simplified version of a real infant gut microbiome study, which tracked the gut microbiota of a newborn over the course of several days after birth. Please have a glance at the paper [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530670/).

### Step 0: Download and inspect the data

Download the tar file with this week's data using the following curl command or by doing a `git pull --no-rebase origin main` in the `~/cmdb-quantbio` directory. If you do the git pull, the tar file should be available in the `~/cmdb-quantbio/assignments/lab/metagenomics/extra_data/` directory. Use `tar -xzvf` to decompress the data file.

```
curl -OL https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/metagenomics/extra_data/metagenomics_data.tgz
tar -xzvf metagenomics_data.tgz
```

Inspect the files within the subfolders of `metagenomics_data` (`step0_givendata`, `step3_givendata`, and `step4_givendata`).

Specifically, within `step0_givendata` you should see

  * a `READS` folder, containing paired-end read data from 8 samples. For simplicity, you can consider these 8 samples to be gut samples collected from the infant from day 0 to day 8 consecutively.

  * You will also see an `assembly.fasta` file, containing the joint assembly of all the samples. Note that entire genomes of the gut bacteria cannot be assembled. Instead, these assembled sequences are the broken up and scrambled pieces of all the various bacterial genomes. You will use this as your reference when you are mapping reads in step 2

  * Finally, there is also a `KRAKEN` folder, containing the approximate classifications of both the reads and the assembled scaffolds. KRAKEN is a very fast and powerful tool that uses k-mer profiles instead of direct alignment to match sequences to the database. Normally you would run KRAKEN on the reads and assembly yourself, but since this required 60GB+ of RAM, they were preprocessed for you. These classifications will be useful to us. Use `less -S` to look over the files and make sure you understand what they mean.

Within `step3_givendata` you should see

  * a `bins` directory that contains 8 fasta files for what we're calling 8 "truth bins" that can be used in step 3 if metaBAT2 doesn't run successfully. They will also be useful in step4 as the given abundance table assumes 8 bins and you might not get 8 bins from metaBAT2. These bin fasta files contain nodes from the assembly that have been grouped together into a bin.

Within `step4_givendata` you should see

  * `abundance_table.tab` are the abundances of each bin in each sample. This table contains the average bin coverage in each sample.

### Step 1: Investigate the taxonomic profile of the reads

We would like to view a detailed taxonomic breakdown of the micro-organisms present in the gut samples for all 8 samples to see how the gut microbiota changes across the days following birth. The KRAKEN classifications of the reads (within `metagenomics_data/step0_givendata/KRAKEN`) should give us this information, but eyeballing the files is not a good idea, so we need a better way to visualize the data.

You should use an existing tool meant for this [Krona-Tools](https://github.com/marbl/Krona/wiki)) to produce pie charts that show the taxonomic breakdown for each sample. (Alternatively, as an advanced exercise, you may get creative with your own code and visualizations). You will need to parse the KRAKEN output to be more usable and then use the [`ktImportText` function in the KronaTools package](https://github.com/marbl/Krona/wiki/Importing-text-and-XML-data) to load and visualize the data.

#### Step 1A: Install Krona-Tools

Use the following code block to install Krona-Tools. You will need to enter your password after using `sudo`.

```
git clone https://github.com/marbl/Krona.git
cd Krona/KronaTools/
sudo ./install.pl
```

#### Step 1B: Parse the KRAKEN output to convert it to a Krona-Tools compatible format

Krona-Tools requires that its input to the `ktImportText` function be tab-delimited separating the taxonomies represented into separate columns. The KRAKEN output files separate the taxonomies with a semicolon rather than a tab. And the first column contains identifier information that we don't need.

Either write and submit your own conversion code or use the following Python code to convert each of the 8 output KRAKEN files:

```
#!/usr/bin/env python

#USAGE: python scriptname.py KRAKEN_sample_input_file.kraken sample_name
#EXAMPLE: python make_krona_compatible.py KRAKEN/SRR492183.kraken SRR492183

import sys

f = open(sys.argv[1])
sample = sys.argv[2]
out = open(sample + "_krona.txt", "w")


for line in f:
    fields = line.strip('\r\n').split(';')
    out.write("\t".join(fields[1:]) + "\n")
```

If you use the provided code, each output file will be named with the sample identifier you provide followed with `_krona.txt`. Record how you called/used this script in a README file.

#### Step 1C: Use Krona-tools

Use the [`ktImportText` function in the KronaTools package](https://github.com/marbl/Krona/wiki/Importing-text-and-XML-data) to load and visualize the converted data. In our conversion, we only provided a list of the taxonomies and not the optional quantities. So be sure to include the `-q` flag.

**Record how you called this function in your README and submit your Krona-Tools html result/visualization.**

**Question 1: In your README, briefly comment on the trends you see in the gut microbiota throughout the first week.**

### Step 2: Deconvolute the assembled scaffolds into individual genomes (binning)

Looking at pie charts is a good way to get a feel for what the microbial communities look like, but we want to know more about individual bacteria! In order to track the changes in abundances of individual prokaryotes throughout time, we need to first extract their genomes from the assembly. As mentioned above, their genomes will be highly fragmented, so the challenge is to "bin" these fragments (contigs or scaffolds) into groups. Ideally, each bin will be a single genome. Our grouping needs to be inclusive enough to get all the parts of each genome, but also specific enough to not include the wrong pieces. Sounds hard, right? It is! Luckily, there is software that can help us...

One software tool that does this is metaBAT2. We've already installed metaBAT2 for you within a conda environment; however, for metaBAT2 to work properly, you will have to do one thing which we will explain at in Step 2C.

Look over the usage examples and guidelines [here for metaBAT2](https://bitbucket.org/berkeleylab/metabat). Note that you need to align the reads from multiple samples to the assembly to get the coverage of each contig in all the samples. This should result in 8 separate alignment files that you will provide to metaBAT2.


**Question 2: In your README, comment on what metrics in the contigs could we use to group them together?**

#### Step 2A: Create an index of the assembly

Use `bwa index` to create an index for the provided assembly (`metagenomics_data/step0_givendata/assembly.fasta`).

#### Step 2B: Use bwa mem to map the reads from each sample to the assembly and samtools sort to create sorted bam output files

Use `bwa mem` to map the reads from each sample to the assembly. We recommend using 4 threads (`-t` argument). Be sure to provide it with the `assembly.fasta` reference file. Recall that the reads for each sample are in the `metagenomics_data/step0_givendata/READS/` directory and that they are paired-end reads so you will need to provide a `_1.fastq` and and `_2.fastq` or each sample.

Then use `samtools` sort on the output from `bwa mem` to output a sorted BAM file that can be given to metaBAT2 as input. Note that you will have 8 separate alignment (sorted BAM) files after Step2B.

#### Step 2C: Repair the metaBAT2 environment

Use this command to repair the metaBAT2 conda environment so that metaBAT2 will successfully run:

```
cp ~/miniconda3/envs/metabat2/lib/libdeflate.so ~/miniconda3/envs/metabat2/lib/libdeflate.0.dylib
```

Activate your metaBAT2 conda environment using the following command.

```
conda activate metabat2
```

##### Step 2D: Run metaBAT2

Because of the way that the metaBAT2 environment is set up, you will have to use the "slightly less easy way" described in the "MetaBAT 2 USAGE: running on command line" section [here in the metaBAT2 documentation](https://bitbucket.org/berkeleylab/metabat/src/master/).

Use the `jgi_summarize_bam_contig_depths` function followed by the `metabat2` function to separate the scaffolds into bins.

**Question 3A: In your README, answer: How many bins did you get?**

**Question 3B: In your README, comment on roughly what percentage of the assembly do they represent?**

(Hint: You can see how many and which "Nodes" are in each bin fasta file by grepping for a `>`, which starts each identification line. Use a pipe and some downstream unix command or python parsing for counting, etc.)

**Question 3C: In your README, comment on whether you think the sizes of each bin look about right, based on what you know about the size of prokaryotic genomes?**

**Question 3D:In your README, describe how you might estimate how complete and how contaminated each bin is?**

### Step 3: Estimate the taxonomy of your putative genomes

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

$ grep NODE_14_length_235766_cov_39.967778 metagenomics_data/step0_givendata/KRAKEN/assembly.kraken

NODE_14_length_235766_cov_39.967778	root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus haemolyticus;Staphylococcus haemolyticus JCSC1435
```

But do this programmatically by tallying up the taxonomic composition of your sample at various levels of taxonomic hierarchy. (The levels of hierarchy are separated by semicolons in the original KRAKEN output)

**Be sure to submit the code/process that you use.**

**Question 4: <br/>(A) Within your README, record your predictions for each bin?<br/>(B) This approach to classification is fast, but not very quantitative. Within your README, propose one method to more robustly infer the taxonomy of a metagenomic bin.**

### Step 4 (Optional): Make a heatmap of the individual bin abundances over time

You have your genomes - now we want to see how their abundances shifted in the baby's gut during these 8 days. Within `metagenomics_data/step4_givendata/` you have an `abundance_table.tab` which contains the abundances of each bin in each sample. These bins correspond to the 8 "truth bins" that an instructor produced. (The truth bins fasta files were provided in `metagenomics_data/step3_givendata/bins/`). This `abundance_table.tab` table contains the average bin coverage in each sample.

Use your favorite software package to make a heatmap. All the samples have the same total read count, so you do not need to standardize your abundance values. Repeat your Step 3 analysis in order to identify the species/genus name that you surmise or predict each bin is best represented by. Then label the bins with these corresponding species/genus name . If you are unsure how to proceed, consider using the seaborn `clustermap` function.

**Submit the heatmap!**

**Question 5: Compare the results of the individual genome abundance analysis to the conclusions you derived from looking at the read taxonomy distributions (from Step 1). Do they agree with each other? What is different?**

## Submit

* All Python code
* Command line record
* Krona Pie chart visualizations from Step 1
* Thoughtful answers to Questions 1-4 in your README
* Optional: Any additional visualizations you made in Step 1 (plus code to produce them)
* Optional: Heatmap from Step 4
* Optional: Thoughtful answers to Question 5 in your README
