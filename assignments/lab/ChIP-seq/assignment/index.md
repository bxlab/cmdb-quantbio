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

The bam files contain reads mapped against the mouse genome mm10 and come from two ChIP-seq replicates targeting Sox2 and their corresponding controls. The BED file is the filtered peak calls from the Klf4 ChIP-seq data. The three bedgraph files are normalized pileups of Klf4 and two H3K27ac ChIP-seq treatment conditions. The D0 and D2 refer to day 0 and day 2 post induction. These datasets include only data for chromosome 17.

#### Mapping reads

Reads have been mapped to the mm10 genome for you for time's sake.

#### Filtering reads

You will be following the protocol described in the paper to filter aligned reads, keeping only those with a quality score of 10 or greater.

- **Filter reads using `samtools view` to only include quality scores >= 10.**

#### Calling peaks

You will use `macs2` to call peaks. `macs2` is preloaded on your computer and has several modes of operation. You will use is called `callpeak`.

Remember that you do not need to provide all of the optional parameters. For `macs2`, you definitely want to provide the target and control (input) samples. It is also typically important to provide a parameter for the effective genome size (since you are only using data from chr17, what is your effective genome size?). Otherwise, the defaults should work reasonably well here. If the enrichment by the antibody was poor, you would need to consider the `mfold` and `extsize` parameters (but more likely, back to the lab). You should also use the `-B` parameter to create a bedgraph of read pileups for the treatment and control conditions. This will be used for plotting in a later step.

`macs2` will write several output files to a directory that you specify. The file containing the peaks has the extension `narrowPeak`. This is an extended BED format, where the first six fields are standard BED. There will also be a pair of files with the read pileups suffixed with `_treat_pileup.bdg` and `control_lambda.bdg`. These include a chromosome, start and stop position, and the normalized number of reads covering the bases in that range (each range is associated with a single score)

- **Run `macs2` to produce a list of peaks for each condition in [BED format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html).**

#### Intersecting peaks

Because the researchers suspected there were many false positives peaks called in their ChIP-seq data, the paper protocol intersected the two replicate peak call sets and kept only those which overlapped. These are the `.narrowPeak` files in your `macs2` output folder. Use `bedtools` to do this. Because you want the intersection or the peaks found in both files, the order of your files doesn’t matter.

- **Intersect peaks and create a file of peaks appearing in both Sox2 replicates.**

#### Colocalization of Sox2 and Klf4

The paper shows that the vast majority of Klf4 peaks are also bound by Sox2. While your numbers will be different since you are only looking at chromosome 17, see if your data give a similar percentage of overlap. Make sure you are using the intersected set of Sox2 peaks. You can use `bedtools` for this task.

- **Find the number of total peaks and overlapping peaks for Klf4 and Sox2 in your data. What is the percentage of Klf4 peaks colocalized with Sox2?**

#### Plot

Next, you will recreate the read pileup tracks from figure 6K. You have been provided with two python scripts in the `assignments/lab/ChIP-seq/extra_data` folder to help with this. One script is meant for scaling the bedgraph files. It does this by adjusting each value such that the mean signal is one read per base. The other script has a function that loads in and bins the bedgraph file data, returning the coordinates of each bin’s midpoint and the sum of read counts for that bin. Further, the returned bin positions and values are cropped to include only those falling in the genomic window shown in the paper’s figure. (Note: when you look at these positions, they will not be the same as what you see in the figure because the paper mapped the reads to the mm9 genome, while you mapped to the more recent build mm10. I have already lifted the coordinates over for you so that the coordinates you see correspond to the paper’s mm9 locations).

While you used both of your replicates to determine true peaks, you will only need one of the treatment bedgraph files for plotting. Combining them in a sensible way is not a productive use of time for this assignment and makes little difference for the region being plotted. You may select the bedgraph from either Sox2 replicate.

Before plotting or cropping the data, you will need to scale the bedgraph files so they are directly comparable with each other. To run the scaling script, use the following command:

```
python scale_bdg.py <original_bdg> <scaled_bdg>
```

Loading the bedgraph files into Python for plotting takes a bit of time, so I’ve also provided you with an `awk` command to crop your data files. Specifically this is cropping/subsetting your data to only include values falling within the plotting window. Since the data are already scaled, you can discard the points outside the plotting range without causing issues.

Everything you’ve done so far focuses on processing the Sox2 data. Now you’ll use your processed Sox2 data as well as the Klf4, day 0 H3K27ac, and day 2 H3K27ac data that I preprocessed for you. To plot, you will want to scale the bedgraph files so they are comparable to each other, crop the scaled files so they don’t take forever to load, load the scaled and cropped data with the provided function, and then produce a 4 panel plot displaying the data.

1. **Scale each of your 4 bedgraph files so they are directly comparable to each other. To run the provided scaling script on one bedgraph file, use the following command:**
  ```
  python scale_bdg.py <input_bdg_filename> <output_scaled_bdg_filename>
  ```
2. **Crop the files using this command for each file:**
  ```
  awk '{ if ($2 < 35502055 && $3 > 35507055) print $0 }' <input_scaled_bdg_filename> > <output_scaled_and_cropped_bdg_filename>
  ```
3. **Create a python script which uses the provided function in bdg_loader.py to load the scaled and cropped bedgraph files**
4. **Within that python script, write code to produce and save a 4 panel plot like the one in figure 6K. Add appropriate track labels. If you want to have the tracks filled in like the original figure, you can look up how to use the *fill_between* function in matplotlib.**

### Part 2: Motif discovery

You will now attempt to discover the sequences that are being bound by under the peaks for Sox2. For this we will use `meme-chip`, which is a suite of many programs. Again, this should already be installed on your computer.

#### Data

You should use your combined Sox2 `.narrowPeak` file for this section.

Also, you want to be able to compare any discovered motifs against known motifs. Download the latest version of [MEME motif databases](https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.23.tgz).

#### Motif finding

You will now find the sequences that are being bound by Sox2. For this we will use `meme-chip`, which is a suite of many programs. Again, this should already be installed on your computer, although it is in its own conda environment. To activate it, use the command:

```
conda activate meme
```

There is a bug in conda with the current version of one of the dependencies so you will also need to run the following command once you have activated the meme environment:

```
conda install -c conda-forge openmpi=4.1.4 -y
```

When you are finished with this assignment, don't forget to exit the environment with `conda deactivate`.

Motif finding is computationally intensive, so you want to enrich for signal. I would suggest passing only the 300 strongest ChIP-seq peaks to MEME (you can pass more, but it will take more time, and passing weak peaks will make the motif(s) more difficult to find). MEME requires the FASTA input for these peaks. We’ve provided the mouse genome FASTA file already on your computer at `/Users/cmdb/data/genomes/mm10.fa`. However, you should not have write access to that directory. Because the file is large and you don't need to alter it, you will be making an alias (or symlink) to it rather than a copy. To do this, use the command:

```
ln -s /Users/cmdb/data/genomes/mm10.fa ./
```

You should now see `mm10.fa` in your current directory.

You’ll use `samtools faidx` to extract the sequences of the signal enriched peaks, however, you’ll need an `awk` command to convert the bed file into the proper format for `samtools faidx`. The `awk` command will take the chromosome name and coordinate range and convert it into a format recognized by `samtools faidx` (as well as the UCSC genome browser).

1. **Sort your intersected Sox2 replicate narrowPeak file by the score column (5th)** 
2. **Keep only the first 300 lines**
3. **Use this command `awk '{ printf "%s:%i-%i\n", $1, $2, $3 }'` to reformat the 300 lines so they can be used as input to `samtools faidx`**
4. **Use `samtools faidx` to extract the sequences of these peaks from the mm10 reference genome. Consider the -r flag for passing in your reformated input**
5. **Run `meme-chip` to perform motif finding in these strongest 300 peaks from Sox2. Consider motif widths up to 7bp (-maxw).**

### Motif identification

Finally, you will compare your discovered motifs against known motifs to see if you can identify transcription factors known to bind to these motifs. You will compare all motifs discovered by MEME to one of the databases you downloaded and then examine your results. The tool for this is `tomtom` which is part of the MEME suite.

1. **You downloaded a large number of databases from the MEME site. Unpack the database file and take a look at what’s inside.**
2. **Using the motif_databases/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme database and your MEME discovered motifs (which can be found in the output folder from `meme-chip` in the file combined.meme), scan your motifs against the database using `tomtom`**
3. **Examine the `tomtom` .html results file in your web browser to see what sorts of matches were found**
4. **Pull out all matches to “KLF4” and “SOX2' from the tomtom.tsv file in the `tomtom` output folder, saving these matches to a separate file.**

#### Submit

* Record of command line commands for the assignment
* File with the number of peaks in Klf4, Sox2, and intersecting peaks
* Script for plotting the track figure
* A pdf of the track figure
* The match profiles from `tomtom` for "KLF4" and "SOX2"


### Advanced exercises

Using the motif position data (accessed through the `meme` html page by clicking on the "Motif Sites in GFF3" link), recreate figure 6B from the paper showing the spacing distribution between Sox2 and Klf4 sites.
