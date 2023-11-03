# Assignment 7: Nanopore Sequencing and Methylation
Assignment Date: Friday, Nov. 5, 2023 <br>
Due Date: Friday, Nov. 12, 2023 @ 10:00am ET <br>

## Assignment Overview

For this lab, you will be working with two datasets. The first is methylation calls across chromosome 2 in the cell line GM24385 using both traditional bisulfite sequencing (using the software Bismark) and direct detection using nanopore sequencing. The second dataset is looking at paired samples for normal and colorectal cancer methylation cells. There are methylation calls for bisulfite sequencing and nanopore sequencing for chr2. You also have the nanopore reads for two subregions of chr2 in the normal and tumor samples.

 The aim of this assignment is to compare methylation calling from the two different technologies and observe the role of tumorigenesis in hypermethylation and allele-specific methylation.

### Part 1: Getting the data

To get the data for this assignment, you will need to download it from [here](https://www.dropbox.com/scl/fi/qz5m5x5v5sr3og3vy9921/ONT_data.tar.gz?rlkey=1bp3j45sqavysnqg2d5n0jxhv&dl=0) (There is a download button in the upper right corner, a down arrow over a horizontal line).

Create a folder for this assignment and move the resulting file `ONT_data.tar.gz` into your newly created folder. Then unpack the file with the following command:

```bash
tar -xzf ONT_data.tar.gz
```

You should have a total of 9 files:

- ONT.cpg.chr2.bedgraph
- bisulfite.cpg.chr2.bedgraph
- normal.ONT.chr2.bedgraph
- normal.ONT.chr2.bam
- tumor.ONT.chr2.bedgraph
- tumor.ONT.chr2.bam
- normal.bisulfite.chr2.bedgraph
- tumor.bisulfite.chr2.bedgraph
- transposons.chr2.bed

### Part 2: Exploring CpG methylation

First, start IGV from the command line (the command is `igv`). Make sure you have hg38 selected as the genome. Then under the file drop-down menu select "Load from File" and add the tracks ONT.cpg.chr2.bedgraph and bisulfite.cpg.chr2.bedgraph. Each of these tracks contain one CpG location per line, the percent of reads with that site methylated, and the read coverage for that site. Explore at different resolutions, comparing the two tracks. Do they match well? Are there regions missing from one track but present in the other?

Q1: Are the majority of the CpG dinucleotides methylated or unmethylated?

### Part 3: Comparing nanopore vs. bisulfite sequencing methylation calling

You will be comparing data from the two methylation calling approaches. This includes answering questions in your `README.md` file as well as creating several plots. Please put all of the plots into a single multi-plot figure with properly labeled axes and titles.

#### Part 3a

Using the above bedgraph files, write a Python script to perform a number of comparisons between the two sets of methylation calls. Your script should be able to:

1. Parse the bedgraph files
2. Calculate the number of sites present only in the bismark file, present only in the nanopore file, and the shared sites as a percentage of total sites (both unique and shared sites) and record them in your `README.md` file.

#### Part 3b:

Plot the distribution of coverages across CpG sites for each track on the same plot. Make sure to indicate which distribution corresponds to which track. In order to visualize both distributions on the same plot, it may be useful to use the `alpha` option to set the transparency of the plotted data.

Q2: How does using nanopore for methylation calling differ from bisulfite sequencing in terms of coverage? Which method appears better and why?

#### Part 3c:

For CpG sites occurring in both bedgraph files, plot the relationship between methylation scores for the two approaches. Because of the number of data points, it is impractical to do this using a scatterplot. Instead, use the numpy function `histogram2d` and plot using the matplotlib function `imshow`. I recommend 100 bins per axis for the histogram as this will make the axis labels match the percent methylation. 

To do this, first run `numpy.histogram2d()` on your data. `numpy.histogram2d()` returns three items - a histogram, the x edges and the y edges. Store the histogram as a variable. Because points are highly concentrated in the corners, it is difficult to see much of the data. Therefore you should transform the histogram using a `log10(data + 1)` transformation. Plot your transformed histogram `imshow()`. 

You also should calculate the Pearson R coefficient for the two sets of methylation calls (non-transformed data). This is easy to do using the numpy function `corrcoef`. Include this value in the title (no more than 3 decimal places, please).

#### Part 3d:

Now, let's examine the matched normal-tumor samples. You should be able to load these bedgraph files with the same parser as the previous files. For each pair of samples (normal and tumor), for all CpG sites in common find the change in methylation (tumor - normal), excluding values with no change. Create a violin plot showing the distribution of methylation changes, one distribution for nanopore and one for bisulfite results. Using common sites between the two approaches, find the Pearson R coefficient for methylation changes and add this value to the title of the plot.

Q3: What can you infer about the two different approaches and their ability to detect methylation changes?
Q4: What is the effect of tumorigenesis on global methylation patterns?

### Part 4: Using IGV to explore differential methylation

Finally, you will be visualizing your data in `IGV`. Before starting, you will need to index the bam files using `samtools`. Load the two bam files and the bismark normal and tumor bedgraph files into the browser. For the bam files, use a two-finger click on the track to pull up a menu, select `Color alignments by` and then select `base modification (5mc)`. This will use the methylation data stored in each read to color methylated sites red and unmethylated sites blue. I recommend also using the same approach to select "squished" to fit the data more easily in the browser.

#### Part 4a:

In the navigation field, type in the gene name "DNMT3A". This should take you to the gene "DNA (cytosine-5)-methyltransferase 3A", a gene responsible for de novo methylation of cytosines and one of the 127 frequently mutated genes identified by the Cancer Genome Atlas project. Find a region of interest and save an image.

Q5: What changes can you observe between the normal and tumor methylation landscape? What do you think the possible effects are of the changes you observed?

#### Part 4b:

In the navigation field, type in the gene name "ZDBF2". This should take you to one of the 91 known imprinted genes. Adjust you field of view to focus only on the first exon of the gene. Do you see much difference between the tumor and normal sample? In order to see the effects of imprinting, you will first need to phase the reads. You can do this by pulling up the menu for each track with a two-finger click and selecting "Cluster (phase) alignments" (accept 2 for the number of clusters).

Q6: What does it mean for a gene to be "imprinted"?
Q7: What is happening when you select the option to phase the reads? What is required in order to phase the reads?

Now that reads are phased, can you identify any allele-specific methylation sites? How does tumorigenesis affect these sites? Is the effect consistent across sites? Save a picture of this region. Drag the view region left or right. Try phasing the reads at different locations. Does it always work?

Q8: Can any set of reads be phased? Explain your answer.

## Submission

For this assignment you should turn in the following:

* Part 3a: Python code to calculate overlap between Bismark and Nanopore (**1pt**)
* Part 3b: Python code to plot distributions of coverages (**1pt**)
* Part 3c: Python code to plot relationship between two approaches and caluclate correlation (**1pt**)
* Part 3d: Python code to compare normal and tumor samples (**1pt**)
* Part 3: Screenshot of multi-plot figure of coverage distributions (3b), heatmap (3c), violin plot (3d) (**2pt**)
* Part 4: 2 images from IGV, one from gene DNMT3A and one from gene ZDBF2 (**1pt**)
* README.md: Answers to questions 1-8: (**3pt**)

## Advanced

One important role of DNA methylation is to inactivate transposons, preventing their spread through the genome and possible disruption of normal function. Included with the data you downloaded is a bed file containing a set of identified transposons on chromosome 2 for the human genome build hg38. Using the ONT.cpg.chr2.bedgraph file, create an average profile of methylation across all of the transposons in the bed file. Include 5kb upstream and downstream of each transposon, binning data into 100bp bins for the flanking regions and dividing each transposon into 100 equal sized bins (thus transposons of different length will have different bin sizes but all have a total of 100 bins covering them). Don't forget to divide the data sums by the binsize to ensure your units are % methylation/bp. Divide the data into quartiles based on the score in the bed file (you will need np.float64 to handle the precision of the scores). You should have one line per quartile of data. Don't forget to include a legend.

The scores in the bed file represent the significance of the sequence similarity to transposon family that each sequence derived from. Hence low scores are younger transposons as few mutations have accumulated. Given that, what does your plot tell you about the relationship of methylation to transposon age?
