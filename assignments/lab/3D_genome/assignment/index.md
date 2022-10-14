# Assignment 6: The 3D Genome
Assignment Date: Friday, Oct. 14, 2022 <br>
Due Date: Friday, Oct. 21, 2022 @ 1:00pm ET <br>

## Lecture

[Lecture slides](https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/3D_genome/slides_asynchronous_or_livecoding_resources/3D_Genome.pdf)

## Assignment Overview

For this lab, you are going to be working with capture Hi-C data from [Nonlinear control of transcription through
enhancer–promoter interactions](https://pubmed.ncbi.nlm.nih.gov/35418676/). In this study, the authors used a transgene in an empty topological domain to explore the effects of enhancer placement both within and outside the domain boundaries.

In today's assignment you will attempt to recreate some of the paper's findings and figures.

### Part 0: Setting up your conda environment

In order to use the HiC-Pro software for Hi-C analysis, you will need to set up a conda environment with all the requirements. This can take some time, so you will get to do this at the **beginning** of class. To do this, create your working folder and copy the file `/Users/cmdb/cmdb-quantbio/assignments/lab/3D_genome/extra_data/environment.yml` in this folder.

Next, you will use this file, which contains a list of programs and versions to install, to create a new conda environment. To do this, use the following command:

```bash
conda env create -f environment.yml -n hicpro
```

### Part 1: Capture Hi-C analysis

#### Data

You will be using capture Hi-C data for two edited cell lines, one with a single deleted CTCF site within the domain of interest and another with two deleted CTCF sites within the domain of interest. Due to time constraints, you will be analyzing a subsampled set of reads. You will also be provided with analyzed data from the complete datasets. You can download the data as follows:

```bash
curl https://bx.bio.jhu.edu/data/msauria/cmdb-lab/3dgenome_data.tar.gz --output 3dgenome_data.tar.gz
curl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes --output mm10.chrom.sizes
```

Next, unpack the datasets.

```bash
tar xzf 3dgenome_data.tar.gz
```

You should now have 13 files:

- mm10.chrom.sizes
- make_plots.sh
- fastq/
	- dCTCF/
    	- SRR14256290_1.fastq
        - SRR14256290_2.fastq
    - ddCTCF/
        - SRR14256299_1.fastq
        - SRR14256299_2.fastq
- mm10.DPNII.frag.bed
- target.bed
- config_hicpro.txt
- matrix/
    - dCTCF_full.6400.matrix
    - ddCTCF_full.6400.matrix
    - 6400_bins.bed
    - dCTCF_full.40000.matrix
    - 40000_bins.bed

The script `make_plots.sh` is a replacement for a buggy one in HiC-Pro.

The fastq folder contains Hi-C sequenced reads organized by sample. The name `dCTCF` corresponds to a single deleted CTCF site, while the name `ddCTCF` corresponds to the CTCF double site deletion.

The file `mm10.DPNII.frag.bed` contains a list of all of the expected genomic fragments produced by cutting the reference genome with the restriction enzyme DpnII, the one used in these experiments. The file `target.bed` contains the region that was enriched for (the capture part of capture Hi-C). The file `config_hicpro.txt` contains settings for running your analysis.

The two .bed files in the matrix folder contain the boundaries of each bin produced by breaking the genome into bins of the size indicated in the file name. The .matrix files contain 3 columns corresponding to the bin # of the first part of an interaction, the bin # of the second part of an interaction, and a normalized score associated with that interaction.

#### Setting up HiC-Pro

The HiC-Pro software is a complete pipeline for analyzing Hi-C data. It performs the following steps:

1. Mapping reads
2. Assigning reads to restriction fragments
3. Quality control
4. Binning reads into uniform-sized genomic bins
5. Normalize binned signal

To run HiC-Pro, you will first need to build it. Start by activating the conda environment you created at the beginning of the assignment:

```bash
conda activate hicpro
```

Next, you will need to clone the HiC-Pro repo. Note that you will be cloning a specific version of the repo rather than the most recent one.

```bash
git clone -b v3.1.0 https://github.com/nservant/HiC-Pro.git
```

Finally, you will need to build HiC-Pro. To do this, move into the repo you just cloned and use the following commands:

```bash
make configure PREFIX=${PWD}/../
make install
```

 This will install HiC-Pro into your working directory. Now exit that directory. Finally, you will need to copy the script `make_plots.sh` into the directory `HiC-Pro_3.1.0/scripts/`

#### Processing the capture Hi-C data

Before running HiC-Pro, you will need to edit the configuration file. There are three lines that you will need to replace `<YOUR_DIRECTORY>` with the complete path to the directory you are working in. You can always see the full path with the command `pwd`. Once you have filled those three entries in, save the config file.

In order to run HiC-Pro, you need to give it the folder where the fastq files are, organized by sample (-i), the name of an output directory to store the results in (-o), and a configuration file (-c). The configuration file was obtained from the repo describing the analysis from the paper and has been updated with the correct paths for various files. You will need to run the program by calling it from `HiC-Pro_3.1.0/bin` which was created in your working directory when you installed HiC-Pro.

Running the analysis will take several minutes. You will see information about each step as it is performed. Once the analysis has finished running, look around in the output directory. You should have five folders. The `bowtie_results` directory has all of the mapped read data. The `hic_results` has all of the data once reads have been assigned to restriction fragments. In this directory, there is a directory `pic` containing several QC plots. Take a look at the `HiCContactRanges` and `HiCFragment` plots (the mapping plots don't tell you anything since only mapped data were selected).

- **What percentage of reads are valid interactions (duplicates do not count as valid)?**
- **What constitutes the majority of invalid 3C pairs? What does it actually mean (you may need to dig into the [HiC-Pro manual](https://github.com/nservant/HiC-Pro/blob/v3.1.0/doc/MANUAL.md))?**


### Part 2: Exploring the heatmaps

#### Creating differential interaction plots

You will now use the data you analyzed as well as that provided at 6400bp resolution to recreate figure 4a from the paper (minus the scale bars). You will have two versions, one with your subsetted data and another with the full data provided. To do this, you will need to take the sparse format that is provided and convert it into a complete matrix for plotting. You have been provided with a starting script for loading the data you will need for this. Your goal is to produce a horizontal 3-panel plot with a heatmap for ddCTCF, dCTCF, and dCTCF - ddCTCF (in this order), covering the region chr15:11170245-12070245 (note that this is different from the paper because they used mm9 and you are using mm10). The heatmaps should come from the `iced` folder and you can use either 6400bp bin file from the `raw` folder in `hic_results/matrix/XXXX`.

To create the plot, you will need to do the following:

1. Filter out data with one or both interaction ends falling outside the desired bin range
2. Log-transform the scores (the dynamic range of data makes it hard to visualize the non-transformed data). Also, shift the data by subtracting the minimum value so the new minimum value is zero (this will prevent issues where there is missing information)
2. Convert the sparse data into a square matrix (note that the sparse data only contains one entry per interaction with the lower-numbered bin in the first column). For one line of the sparse data format, the data relates to the full matrix as follows:

	```python
	mat[sparse['F1'][i], sparse['F2'][i]] = sparse['score'][i]
	```
	
4. Plot the two matrices using the same maximum value (set vmax in `imshow`). I suggest using the `magma` color map, although you need to flip your scores to mimic the paper figure
5. For the difference plot, I suggest using the `seismic` color map and `norm=colors.CenteredNorm`. It helps to remove the distance dependent signal and smooth the data first as there is noise. You can use the following function to remove the distance dependent signal:

	```python
	def remove_dd_bg(mat):
	    N = mat.shape[0]
	    mat2 = numpy.copy(mat)
	    for i in range(N):
	        bg = numpy.mean(mat[numpy.arange(i, N), numpy.arange(N - i)])
	        mat2[numpy.arange(i, N), numpy.arange(N - i)] -= bg
	        if i > 0:
	            mat2[numpy.arange(N - i), numpy.arange(i, N)] -= bg
	    return mat2
	 ```

 	You can use this function to create smoothed matrices before subtracting:

	```python
	def smooth_matrix(mat):
	    N = mat.shape[0]
	    invalid = numpy.where(mat[1:-1, 1:-1] == 0)
	    nmat = numpy.zeros((N - 2, N - 2), float)
	    for i in range(3):
	        for j in range(3):
	            nmat += mat[i:(N - 2 + i), j:(N - 2 + j)]
	    nmat /= 9
	    nmat[invalid] = 0
	    return nmat
	```

- **Were you able to see the highlighted difference from the original figure?**
- **What impact did sequencing depth have?**
- **What does the highlighted signal indicate?**

#### Finding insulation scores

Next you will need to score the provided 40kb resolution data to determine the level of insulation between each bin. You will use the whole range of data for this region which spans bins 54878 to 54951. To do this, you will first need to log-transform the data as before, subtracting the minimum score after transformation. You will then need to convert the data into a matrix.

Once you have a matrix, you will need to find the insulation score by taking the mean of the 5x5 square of interactions between the 5 upstream bins and 5 downstream binsaround the target.

You can use the following syntax to find the mean of a 5x5 block. This one is for upstream of the test point *i*.

```python
numpy.mean(mat[(i - 5):i, i:(i + 5)])
```

Note that the insulation score is found at the boundary between two bins (if i==5, then the insulation score corresponds to the start of bin 5, assuming that bin 0 is the first bin).

#### Plotting insulation scores

Finally, plot the heatmap for bins 54878 to 54951 (chr15:10400000-13400000). Use the log-transformed data but not distance-corrected. You should also plot the insulation scores below the heatmap, lining them up with the heatmap. I found that using the following worked nicely for this purpose:

```python
fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(5,6.25))
ax[0].axis('off')
plt.margins(x=0)
ax[1].set_xlim(10400000, 13400000)
plt.subplots_adjust(left=0.15,
                bottom=0.1,
                right=1.0,
                top=1.0,
                wspace=0.4,
                hspace=0.0)
```

## Submission

For this assignment you should submit four things:
1. A README.md with answers to the questions
2. A three panel plot for your analyzed data
3. A three panel plot for the full data
4. A two panel plot with insulation scores


## Advanced exercise

Find peaks in your insulation scores and add lines to the heatmap corresponding to these peaks. You can define peaks as points with values higher than both neighboring points.
