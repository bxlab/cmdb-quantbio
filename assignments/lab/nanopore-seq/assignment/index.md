# Assignment 7: Nanopore Sequencing and Methylation
Assignment Date: Friday, Oct. 28, 2022 <br>
Due Date: Friday, Nov. 4, 2022 @ 1:00pm ET <br>

## Lecture

[Lecture slides](https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/nanopore-seq/slides_asynchronous_or_livecoding_resources/Nanopore.pdf)

## Assignment Overview

For this lab, you will be working with nanopore data from Oxford Nanopore. The data were basecalled to identify methylated cytosine in addition to standard sequence. Reads for four regions have been extracted along with their corresponding methylation calls. If you look at the entries in the bam file, you will see that there is an additional pair of tags, MM and ML, followed by a series of comma-separated numbers. These values base counts between methylated bases and call probabilities, respectively.

 The aim of this assignment is to identify differential methylation at several known imprinted genes from the paper [The role of imprinted genes in humans](https://pubmed.ncbi.nlm.nih.gov/22771538/). If you observe figure 2, you can see a series of regions with differential methylation.

 ![Figure 2](./fig2.png)

For each of these regions, create an image from IGV showing reads separated by parent and highlighting their methylation patterns. There are several steps needed to achieve this.

1. Call and phase variants for each region
2. Mark reads with the correct haplotype tag
3. Split reads into two files base on their haplotype
4. Load data into IGV and set up methylation-specific viewing parameters
5. Locate differentially-methylated feature(s) for each region
6. Save an image of each region

### Part 1: Setting up your conda environment

In order to call variants and phase them, you will be using two programs, Medaka and whatshap. You will need to create a conda environment to do this.

```bash
conda env create -n Medaka medaka pandas pysam gradle openjdk=11
```

This will install most of what you need. Once you have activated this environment, you will need to replace the current version of whatshap, as there is an incompatibility with the installed version. Unfortunately you can't do this with `conda` so you will need to use the following command instead:

```bash
pip install whatshap==1.0
```

This will force an updated version of whatshap to be installed. When you run `medaka`, it may complain that whatshap is not the correct version but will run anyway. Ignore these warnings.

You will also need to install a new version of `IGV` for this assignment as the functionality for viewing methylation data without reprocessing bam files was only recently added. To do this, you will be cloning the `IGV` repo and building the application.

```bash
git clone https://github.com/igvteam/igv.git
```

Once you have cloned the repo, change into the igv directory and use the following command to build the program:

```bash
cd igv
./gradlew createDist
cd ../
```

This will create a new executable in the folder `igv/build/IGV-dist/igv.sh`. For convenience, I suggest creating a symlink (like an alias) in your homework directory.

```bash
ln -s ${PWD}/igvbuild/IGV-dist/igv.sh ./
```

Now, when you want to start `IGV`, simply type `./igv.sh`.

### Part 2: Getting the data

All of the data is downloadable in a compressed archive. To download and unpack it, use the command:

```bash
curl https://bx.bio.jhu.edu/data/msauria/cmdb-lab/ont_data.tar.gz --output ont_data.tar.gz
tar xzf ont_data.tar.gz
```

This will give you 5 files:

- methylation.bam - your methylation-called reads
- methylation.bam.bai - the index for your bam file
- methylation.bed.gz - a compressed bed file containing methylation sites
- methylation.bed.gz.tbi - the index for the compressed bed file
- regions.bed - a bed file containing the regions that were used to filter the data

You will also need a copy of the human genome, since using it directly from your data directory will complain about write priviledges.

```bash
cp /Users/cmdb/data/genomes/hg38.fa ./
```

### Part 3: Calling variants and haplotypes

In order to call variants, you will be using `medaka_variant`. In order to do this effectively, medaka needs to know which model was used to do the basecalling, both for variant calling and phasing. You should use `r941_prom_high_g360` for the arguments `-s` and `-m`. You will also need the phased vcf file as output. Only one region can be specified for Medaka at a time so you will need to generate a phased vcf file for each region in `regions.bed` using the format `chr:start-end`.

### Part 4: Phasing the read data

The next step is to assign reads a haplotype and then split the bam file by haplotype. To do this, you will be using `whatshap haplotype` and `whatshap split`, respectively. As in the previous step, you will be doing this for each region independently. There is an error in the `whatshap` usage information about the format of specifying a region. You need to use the format `chr:start:end`. You will also need to output a haplotag list for the next step. Finally, note that the VCF file needed for input should be the compressed and indexed one from the medaka output (i.e. ending in `.gz`)

Once you have tagged the reads, you will need to split the file into two bam files, one for each haplotype. Use `whatshap split` and the haplotag list that you just generated.

### Part 5: Visualization

Finally, you will be visualizing your data in `IGV` (the new one you just installed). While you can load the haplotype-specific bam files region by region, it is easy to concatenate them into one bam file per haplotype using `samtools cat`. Once you have loaded the haplotype-specific bam files, you will need to change the coloring scheme to display the embedded methylation data. To do this, two-finger click on each track and then select `Color alignments by->base modification (5mC)`.

For each region, you should locate one of the differentially-methylated features shown in the above figure. Zoom in and center the feature in the IGV window such that it takes up approximately the middle fifth of the window.

## Submission

For this assignment you should turn in a record of commands you used to do the assignment and one igv image for each region (4 in all).