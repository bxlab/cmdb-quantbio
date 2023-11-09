# The 3D Genome

## Assignment Overview

For this lab, you will be using Promoter-Capture Hi-C (pCHiC) data from the human cell line GM12878. This cell line is derived from B-cells and has been highly studied, meaning that there are many datasets available for it. In this assignment, you will be running the software [Chicago](https://bitbucket.org/chicagoTeam/chicago/src/master/) to analyze interaction counts to find significant promoter interactions. All of the data are mapped to the hg19 build of the human genome.

The aim of this assignment is to be able to call significant interaction peaks in pCHiC data and visualize these interactions using the UCSC genome browser. Then, using additional genome annotations identify potential regulatory elements for relevant genes in this cell line.

## Data

To get the data, download it from [here](https://www.dropbox.com/scl/fi/7d33oefr0633jpf8yufbn/pchic_data.tar.gz?rlkey=663s185pzdolx1xnj9hvndnfk&dl=0) and move it into your submission folder for this assignment. Then unpack it with the command

```bash
tar -xzf pchic_data.tar.gz
```

You should have an R script called `runChicago.R` and a folder called `raw` with three sub-folders, `Design`, `Features`, and `PCHIC_data`. There are several important files that you will need for the assignment.

1. `raw/Design/h19_chr20and21.rmap` contains the genomic coordinates and fragment number (just a label for the fragment) for each restriction digest fragment in your experiment.
2. `raw/Design/h19_chr20and21.baitmap` contains the genomic coordinates, fragment number, and gene name for each promoter bait fragment. This is essentially a subset of the `h19_chr20and21.rmap` file, containing only those restriction digest fragments for which you had a corresponding bait in your pCHiC experiment.
3. `raw/Features/` contains a set of `.narrowPeak` files from previous ChIP-seq experiments in this cell line. There are `.narrowPeak` files for several different interesting regulatory features like CTCF and H3K4me3. You'll be intersecting your pCHiC results with these regulatory features to learn more about how these features regulate long-range interactions. The `raw/Features/featuresGM.txt` file lists the name of each regulatory feature and its corresponding `.narrowPeak` file, and will be needed when running `chicago`.
4. `raw/PCHIC_data/` contains the raw number of reads for each pair of fragments across 3 replicates. This is essentially the raw output of your pCHiC experiment and `chicago` will use this raw data to identify significant interactions.

## Setting up your analysis environment

You will need to use your **x86** version of iTerm for this assignment. Open `iTerm-x86` and run the following command to create and activate your conda environment:

```bash
mamba create -n chicago bioconductor-chicago r-argparser matplotlib
mamba activate chicago
```

You should now see the environment indicator `(chicago)` to the left of your terminal prompt. Remember when you're finished, to get out of this conda environment you can close the window or run the command `mamba deactivate`.

## Exercises

### Exercise 1: Running `Chicago`

The data have already been aligned and assigned to specific restriction fragments. The next step is to determine which fragment pairs are cooccurring at a significantly higher frequency than expected. To do this, you will use the R package `Chicago`. However, rather than running this using the R GUI, you will run this as a script from the command line. You can call the analysis script using the command `Rscript runChicago.R`.

You will need to fill in several options, including the design directory, the feature list file, and each of the sample replicates. Record the command you used in your `README.md` file. Once you've run this you should have a folder with the output prefix you selected containing a variety of output files.

### Exercise 2: Feature enrichment

Inside the output folder there is a subfolder named `enrichment_data`. Examine the feature overlaps (the file without a number is the enrichment across all replicates and the one we're interested in).

Q1. Do these enrichments make sense to you? Are any surprising? Explain your reasoning briefly for each feature.

### Exercise 3: Creating a custom track

`Chicago` outputs a file of significant interactions but in a format specific for the WashU genome browser, which we will not be using. Take a look at the file, which should be `<OUTNAME>/data/<OUTNAME>_washU_text.txt` which `<OUTNAME>` is the output prefix you used when running `Chicago`. The format that you will need to visualizing in the UCSC genome browser can be found [here](https://genome.ucsc.edu/goldenPath/help/interact.html).

Note that the first set of coordinates are the smallest and largest for the pair of interactions (i.e. the start and end will come from opposite ends of the interaction). Use the bait list (`raw/Design/h19_chr20and21.baitmap`) to add gene names to interactions where appropriate (every interaction will have at least one end from a promoter bait fragment, some will have both ends when two promoters are interacting). You should also find the maximum score and scale values from 0-1000 for the `score` field. This will allow you to visualize strong vs. weak interactions in the browser.

You will also want to add the following line at the beginning of your interaction bed file so it properly displays:

`track type=interact name="pCHIC" description="Chromatin interactions" useScore=on maxHeightPixels=200:100:50 visibility=full`

You can load a session in the UCSC genome browser with [this link](https://genome.ucsc.edu/s/msauria/3D_genome_assignment). It has a variety of annotation tracks for GM12878 already loaded. In order to add your custom track, use the button "add custom tracks" under the browser display and then above the "Paste URLs or data" field, use the button "Choose File" to select the interaction bed file you just created.

### Exercise 4: Finding the most significant promoter-promoter and promoter-enhancer interactions

Using either the data from the `washU_text.txt` file or your interaction bed file, Find the 6 top-scoring interactions between two promoter bait fragments and the top 6 interactions between promoter bait and non-bait fragments. Make sure to include the gene names for the bait fragments.

Select at least two of the top 6 promoter-enhancer interactions and navigate to the associated gene in the UCSC genome browser. You may want to click on your custom track in the track selection section and set the minimum score to somewhere between 300-500 to reduce the number of interactions shown. Find likely enhancer-promoter interactions based on the interaction track and annotations. Focus the browser range to only include both ends of these interactions and save an image (if you two-finger click on the browser window, a menu pops up with the option "View image". This will open a new tab which you can two-finger click and select "Save Image As...").

In addition to finding likely enhancers for your target gene, look it up on the [NCBI website](https://www.ncbi.nlm.nih.gov/gene/). You may need to click on the gene link once you search to get all of the information. Specifically, look at the distribution of expression across tissues.

Q2. Does it make sense for this gene to be interacting with enhancers? Explain.

## Submission

You will need to turn in the following:

1. A README.md file with the command you used to run `Chicago` and answers to the two questions.
2. The interaction bed file (git will complain about adding this since it is a bed file which is in your .gitignore list. Read carefully to make sure you actually add it to your commit).
3. Your two lists of most significant interactions, promoter-promoter and promoter-enhancer.
4. Your two UCSC genome browser screen shots of your selected promoter-enhancer interaction(s).

## Advanced Exercise 1: Plotting top interactions (OPTIONAL)

Write a script that takes in the raw data and significant interaction file as well as a target gene name and plots all of the interactions for that gene as a scatterplot of distance from fragment centers between the promoter bait and other fragment by the number of reads. Color the points by the significance, black to red, for lowest to highest, respectively. For interactions not in the significant interaction file, give them a score of zero. You can limit the x-axis range to +/- 1Mb. Plot interaction profiles for at least two genes in your top 6 lists. 

<br><br>
