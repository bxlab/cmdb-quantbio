# Assignment 5: ChIP-seq & Motif Finding
Assignment Date: Friday, Oct. 7, 2022 <br>
Due Date: Friday, Oct. 14, 2022 @ 1:00pm ET <br>

## Lecture

[Lecture slides](https://github.com/bxlab/qbb2021/raw/main/week5/Epigenetics%20-%20ChIPseq.pdf)

## Assignment Overview

For this lab, we are going to use some datasets from [Sox2 and Klf4 as the Functional Core in Pluripotency
Induction without Exogenous Oct4](https://pubmed.ncbi.nlm.nih.gov/31722212/). In this study, the authors demonstrate that contrary to what was commonly believed, *Oct4* is not necessary for the induction of cells back to a pluripotent state. However, *Sox2* and *Klf4* are neccessary and must be in balanced stochiometry in order to produce this induction.

In today's assignment you will attempt to recreate some of the paper's findings and figures.

### Part 1: ChIP-seq analysis

#### Data

You will be using ChIP-seq data for Sox2, Klf4, and H3K27ac. For the sake of time and file size, you will only be processing Sox2 treatment and input mapped reads from chromosome 17. You will also have peak calls appearing in both replicates for day2 Klf4 and read pileups for day2 Klf4 (1st replicate), day0 H3K27ac, and day2 H3K27ac. You can download the data as follows:

```bash
curl https://bx.bio.jhu.edu/data/msauria/cmdb-lab/chipseq_data.tar.gz --output chipseq_data.tar.gz
```

Next, unpack the datasets.

```bash
tar xzf chipseq_data.tar.gz
```

You should now have 8 files:

- D2_Sox2_R1.bam
- D2_Sox2_R1_input.bam
- D2_Sox2_R2.bam
- D2_Sox2_R2_input.bam
- D2_Klf4_peaks.bed
- D2_Klf4_treat.bdg
- D0_H3K27ac_treat.bdg
- D2_H3K27ac_treat.bdg

These datasets include only data for chromosome 17.

#### Mapping reads

Reads have been mapped to the mm10 genome for you for time's sake.

#### Filtering reads

You will be following the protocol described in the paper so filter aligned reads, keeping only those with a quality score of 10 or greater.

1. Filter reads using `samtools view` to only include quality scores >= 10.

#### Calling peaks

You will use `macs2` to call peaks which is preloaded on your computer. `macs2` has several modes of operation, but the one you will use is called `callpeak`.

Remember that you do not need to provide all of the optional parameters. For `macs2`, you definitely want to provide the target and control (input) samples. It is also typically important to provide a parameter for the effective genome size (since you are only using data from chr17, what is your effective genome size?). Otherwise, the defaults should work reasonably well here. If the enrichment by the antibody was poor, you would need to consider the `mfold` and `extsize` parameters (but more likely, back to the lab). You should also use the `-B` parameter to create a bedgraph of read pileups for the treatment and control conditions.

`macs2` will write several output files to a directory that you specify. The file containing the peaks has the extension `narrowPeak`. This is an extended BED format, where the first six fields are standard BED. There will also be a pair of files with the read pileups suffixed with `_treat_pileup.bdg` and `control_lambda.bdg`. These include a chromosome, start and stop position, and the normalized number of reads covering the bases in that range (each range is associated with a single score)

2. Run `macs2` to produce a list of peaks for each condition in [BED format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html).

#### Intersecting peaks

Because there were many false positive peaks in their ChIP-seq data, the paper protocol intersected the two replicate peak calls, keeping only those that overlapped. You should be pretty familiar with doing this with `bedtools` now. It doesn't matter which file you keep peaks from, R1 or R2.

3. Intersect peaks and create a file of peaks appearing in both Sox2 replicates.

#### Colocalization of Sox2 and Klf4

The paper shows that the vast majority of Klf4 peaks are also bound by Sox2. While your numbers will be different since you are only looking at chromosome 17, see if your data give a similar percentage of overlap. Make sure you are using the intersected set of Sox2 peaks.

4. Find the number of total peaks and overlapping peaks for Klf4 and Sox2 in your data. What is the percentage of Klf4 peaks colocalized with Sox2?

#### Plot

Next, you will recreate the plot in figure 6K (the read pileup tracks only). You have been provided with two python scripts in the `assignments/lab/ChIP-seq/extra_data` folder. These are a python script for scaling the bedgraph files for proper display and a python function for loading in and binning the bedgraph files, returning the coordinates of each bin's midpoint and the sum of read counts for that bin. The bin positions and values are cropped to only those falling in the genomic window shown in the paper's figure. Note that the coordinates you will be looking at are not the same as in the figure as the paper mapped data to mm9 while you are using mm10. I have already lifted the coordinates over for you.

Before plotting or cropping the data, you will need to scale the bedgraph files so they are directly comparable with each other. This is done by adjust each value such that the mean signal is one read per base. To run the scaling script, use the following command:

```
python scale_bdg.py <original_bdg> <scaled_bdg>
```

Loading these files take a bit of time, so I suggest you create a cropped version of each of the bedgraphs after scaling but before plotting. You can crop the files using the command:

```
awk '{ if ($2 < 35502055 && $3 > 35507055) print $0}'
```

5. Produce a 4 panel plot like the one in figure 6K with appropriate scaling and track labels.


### Part 2: Motif discovery

You will now attempt to discover the sequences that are being bound by under the peaks for Sox2. For this we will use `meme-chip`, which is a suite of many programs. Again, this should already be installed on your computer.

#### Data

You should use your combined Sox2 `.narrowPeak` file for this section.

Also, we want to be able to compare any discovered motifs against known motifs. Download the latest version of [MEME motif databases](https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.23.tgz).

#### Motif finding

Motif finding is computationally intensive, so you want to enrich for signal. I would suggest passing only the 300 strongest ChIP-seq peaks to MEME (you can pass more, but it will take more time, and passing weak peaks will make the motif(s) more difficult to find). This means you will need to sort by the score column (5th) and keep only the first 300 lines. You will need to get the sequences corresponding to the peaks you have found. The mouse genome FASTA file is already on your computer at `/Users/cmdb/data/genomes/mm10.fa`. You can use `samtools faidx` to extract sequence from a FASTA file but you need a specially formatted region file. You can use the following `awk` command to convert a bed file into the proper format:

```
awk '{ printf "%s:%i-%i\n", $1, $2, $3 }'
```

6. Using `meme-chip`, perform motif finding in the strongest 300 peaks from Sox2. Consider motif widths up to 7bp (-maxw).

### Motif identification

Finally, compare your motifs against known motifs. You downloaded a large number of databases from the `Meme` site. Unpack the database file and take a look at what's inside. You will be specifically using the `motif_databases/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme` database. Your motifs can be found in the output folder from `meme-chip` in the file combined.meme. Pass the whole file (all motifs will be compared to the database). Examine the `tomtom` results file (the .html file) in your web browser to see what sorts of matches were found.

7. Scan your motifs against the above database using `tomtom` (part of the meme suite) to find matches to known motifs.

8. Pull out all matches to "KLF4" and "SOX2" from the tomtom.tsv file in the output folder and save them to a separate file.

#### Submit

* Record of command line commands for the assignment
* File with the number of peaks in Klf4, Sox2, and intersecting peaks
* Script for plotting the track figure
* A pdf of the track figure
* The match profiles from `tomtom` for "KLF4" and "SOX2"


### Advanced exercises

Using the motif position data (accessed through the meme html page by clicking on the "Motif Sites in GFF3" link), recreate figure 6B from the paper showing the spacing distribution between Sox2 and Klf4 sites.
