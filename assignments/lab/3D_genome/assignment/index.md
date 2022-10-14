<header>
<h1> Assignment 6: The 3D Genome</h1>
</header>

Assignment Date: Friday, Oct. 14, 2022 <br>
Due Date: Friday, Oct. 21, 2022 @ 1:00pm ET <br>

<header>
<h2> Lecture</h2>
</header>

<a href="https://github.com/bxlab/cmdb-quantbio/raw/main/assignments/lab/3D_genome/slides_asynchronous_or_livecoding_resources/3D_Genome.pdf">Lecture slides</a>


<header>
<h2> Assignment Overview</h2>
</header>

For this lab, you are going to be working with capture Hi-C data from <a href="https://pubmed.ncbi.nlm.nih.gov/35418676/">Nonlinear control of transcription through
enhancer–promoter interactions</a>. In this study, the authors used a transgene in an empty topological domain to explore the effects of enhancer placement both within and outside the domain boundaries.

In today's assignment you will attempt to recreate some of the paper's findings and figures.

<header>
<h3>Part 1: Getting, Exploring, and Commenting on the HiCPro analysis results</h3>
</header>

<header>
<h4>Get the results</h4>
</header>

Due to issues running HiCPro, please download the analysis results from running HiCPro on the subsampled data:


<pre><code>
curl https://bx.bio.jhu.edu/data/msauria/cmdb-lab/hicpro_analysis.tar.gz --output hicpro_analysis.tar.gz
tar -xzvf hicpro_analysis.tar.gz
</code></pre>

Please also download the analysis results from Mike running HiCPro on the full dataset:

<pre><code>
curl https://bx.bio.jhu.edu/data/msauria/cmdb-lab/3dgenome_data.tar.gz --output 3dgenome_data.tar.gz
tar -xzvf 3dgenome_data.tar.gz
</code></pre>

<header>
<h4> Explore the results </h4>
</header>

After unzipping the <code>hicpro_analysis.tar.gz</code> tar file, you will have a directory called <code>analysis</code> with 5 subfolders

<ol>
	<li> <code>bowtie_results</code></li>
	<li> <code>hic_results</code></li>
	<li> <code>fastq</code></li>
	<li> <code>logs</code></li>
	<li> <code>tmp</code></li>
</ol>

The <code>bowtie_results</code> directory has all of the mapped read data. The <code>hic_results</code> has all of the data once reads have been assigned to restriction fragments. In this directory, there is a directory <code>pic</code> containing several QC plots. Take a look at the <code>HiCContactRanges</code> and <code>HiCFragment</code> plots (the mapping plots don't tell you anything since only mapped data were selected).

You will want to focus on the <code>hic_results</code> directory and its subdirectories <code>pic</code> and <code>matrix</code> for the files you'll use as input when plotting and

After unzipping the <code>3dgenome_data.tar.gz</code> tar file, you will have 15 files.

<ul>
<li><code>make_plots.sh</code></li>
<li><code>load_data.py</code></li>
<li><code>fastq/</code></li>
	<ul>
		<li><code>dCTCF/</code></li>
			<ul>
    		<li><code>SRR14256290_1.fastq</code></li>
        <li><code>SRR14256290_2.fastq</code></li>
			</ul>
    <li><code>ddCTCF/</code></li>
      <ul>
				<li><code>SRR14256299_1.fastq</code></li>
        <li><code>SRR14256299_2.fastq</code></li>
			</ul>
		</ul>
<li><code>mm10.DPNII.frag.bed</code></li>
<li><code>target.bed</code></li>
<li><code>config_hicpro.txt</code></li>
<li><code>environment.yml</code></li>
<li><code>matrix/</code></li>
	<ul>
    <li><code>dCTCF_full.6400.matrix</code></li>
    <li><code>ddCTCF_full.6400.matrix</code></li>
    <li><code>6400_bins.bed</code></li>
    <li><code>dCTCF_full.40000.matrix</code></li>
    <li><code>40000_bins.bed</code></li>
	</ul>
</ul>

You only need to focus on the <code>matrix</code> directory and the python script <code>load_data.py</code>

<header>
<h5><code>matrix</code> directory</h5>
</header>

This directory contains the analysis results from running HiCPro on the full dataset for both the dCTCF and the ddCTCF genotypes. The two <code>.bed</code> files in the <code>matrix</code> folder contain the boundaries of each bin produced by breaking the genome into bins of the size indicated in the file name. The <code>.matrix</code> files contain 3 columns corresponding to the bin # of the first part of an interaction, the bin # of the second part of an interaction, and a normalized score associated with that interaction.

<header>
<h5><code>load_data.py</code></h5>
</header>

You will want to use this script to read the matrix and bin files when you plot the heatmaps in Part 2 of the assignment.

<header>
<h4>Comment on the results</h4>
</header>

Considering the plots in the <code>analysis/hic_results/pic</code> directory, comment on the following:

<ul>
	<li><b>What percentage of reads are valid interactions (duplicates do not count as valid)?</b></li>
	<li><b>What constitutes the majority of invalid 3C pairs? What does it actually mean (you may need to dig into the <a href="https://github.com/nservant/HiC-Pro/blob/v3.1.0/doc/MANUAL.md">HiC-Pro manual</a>)?</b></li>
</ul>

You may find this post helpful: <a href="https://nservant.github.io/HiC-Pro/RESULTS.html">https://nservant.github.io/HiC-Pro/RESULTS.html</a>

<header>
<h3>Part 2: Exploring the results by plotting heatmaps</h3>
</header>


<header>
<h4>Creating differential interaction plots</h4>
</header>

You will now use the subsampled data which was analyzed as well as the full data (both provided at 6400bp resolution) to recreate figure 4a from the paper (minus the scale bars). Specifically you will plot two different versions, one with the subsetted data and another with the full data provided.

Your goal is to produce a horizontal 3-panel plot with a heatmap for ddCTCF, dCTCF, and dCTCF - ddCTCF (in this order), covering the region chr15:11170245-12070245 (note that this is different from the paper because they used mm9 and you are using mm10).

You have been provided with a starting script for loading the data you will need for this. You will input

<ul>
	<li> two sparse format matrices </li>
	<ul>
		<li>one ddCTCF</li>
		<li>one dCTCF</li>
	</ul>
	<li>one bin file</li>
	<li>the output name for your heatmap</li>
</ul>

You will want to input the sparse format that is provided as results from HiCPro. The sparse format data you want should come from the <code>iced</code> or normalized data folder. For the input bin file, you can use either 6400bp bin file from the <code>raw</code> folder in <code>hic_results/matrix/XXXX</code>.

To create the plot, you will need to do the following:

<ol>
	<li>Filter out data with one or both interaction ends falling outside the desired bin range</li>
	<li>Log-transform the scores (the dynamic range of data makes it hard to visualize the non-transformed data).</li>
	<li>Also, shift the data by subtracting the minimum value so the new minimum value is zero (this will prevent issues where there is missing information)</li>
	<li>Convert the sparse data into a square matrix (note that the sparse data only contains one entry per interaction with the lower-numbered bin in the first column). By converting the sparse matrix it into a complete matrix for plotting, you have two entries per interaction. For one line of the sparse data format, the data relates to the full matrix as follows:</li>

	<pre><code>
	mat[sparse['F1'][i], sparse['F2'][i]] = sparse['score'][i]
	</code></pre>

	<li>Plot the two matrices using the same maximum value (set <code>vmax</code> in <code>imshow</code>). I suggest using the <code>magma</code> color map, although you need to flip your scores to mimic the paper figure</li>
	<li>For the difference plot, I suggest using the <code>seismic</code> color map and <code>norm=colors.CenteredNorm</code>. It helps to remove the distance dependent signal and smooth the data first as there is noise. You can use the following function to remove the distance dependent signal:</li>
</ol>

	<pre><code>
	def remove_dd_bg(mat):
	    N = mat.shape[0]
	    mat2 = numpy.copy(mat)
	    for i in range(N):
	        bg = numpy.mean(mat[numpy.arange(i, N), numpy.arange(N - i)])
	        mat2[numpy.arange(i, N), numpy.arange(N - i)] -= bg
	        if i > 0:
	            mat2[numpy.arange(N - i), numpy.arange(i, N)] -= bg
	    return mat2
	 </code></pre>

 	You can use this function to create smoothed matrices before subtracting:

	<pre><code>
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
	</code></pre>

<ol>
	<li><b>Were you able to see the highlighted difference from the original figure?</b></li>
	<li><b>What impact did sequencing depth have?</b></li>
	<li><b>What does the highlighted signal indicate?</b></li>
</ol>

<header>
<h4> Finding insulation scores </h4>
</header>

Next you will need to score the provided 40kb resolution data to determine the level of insulation between each bin. You will use the whole range of data for this region which spans bins 54878 to 54951. To do this, you will first need to log-transform the data as before, subtracting the minimum score after transformation. You will then need to convert the data into a matrix.

Once you have a matrix, you will need to find the insulation score by taking the mean of the 5x5 square of interactions between the 5 upstream bins and 5 downstream bins around the target.

You can use the following syntax to find the mean of a 5x5 block. This one is for upstream of the test point <i>i</i>.

<pre><code>
numpy.mean(mat[(i - 5):i, i:(i + 5)])
</code></pre>

Note that the insulation score is found at the boundary between two bins (if i==5, then the insulation score corresponds to the start of bin 5, assuming that bin 0 is the first bin).

<header>
<h4>Plotting insulation scores </h4>
</header>

Finally, plot the heatmap for bins 54878 to 54951 (chr15:10400000-13400000). Use the log-transformed data but not distance-corrected. You should also plot the insulation scores below the heatmap, lining them up with the heatmap. I found that using the following worked nicely for this purpose:

<pre><code>
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
</code></pre>

<header>
<h2>Submission</h2>
</header>

For this assignment you should submit four things:

<ol>
	<li>A README.md with answers to the questions</li>
	<li>A three panel plot for the subsampled analyzed data</li>
	<li>A three panel plot for the full analyzed data</li>
	<li>A two panel plot with insulation scores</li>
</ol>

<header>
<h2>Advanced exercise</h2>
</header>

Find peaks in your insulation scores and add lines to the heatmap corresponding to these peaks. You can define peaks as points with values higher than both neighboring points.

***
<details><summary><b>What you would have done to run HiCPro:</b></summary>

<header>
<h3>Part 0: Setting up your conda environment</h3>
</header>

In order to use the HiC-Pro software for Hi-C analysis, you will need to set up a conda environment with all the requirements. This can take some time, so you will get to do this at the <b>beginning</b> of class. To do this, create your working folder and copy the file <code>/Users/cmdb/cmdb-quantbio/assignments/lab/3D_genome/extra_data/environment.yml</code> in this folder.

Next, you will use this file, which contains a list of programs and versions to install, to create a new conda environment. To do this, use the following command:

<code>
conda env create -f environment.yml -n hicpro
</code>

<header>
<h3> Part 1: Capture Hi-C analysis </h3>
</header>

<header>
<h4> Data </h4>
</header>

You will be using capture Hi-C data for two edited cell lines, one with a single deleted CTCF site within the domain of interest and another with two deleted CTCF sites within the domain of interest. Due to time constraints, you will be analyzing a subsampled set of reads. You will also be provided with analyzed data from the complete datasets. You can download the data as follows:

<code>
curl https://bx.bio.jhu.edu/data/msauria/cmdb-lab/3dgenome_data.tar.gz --output 3dgenome_data.tar.gz
</code>
<br/><br/>
<code>
curl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes --output mm10.chrom.sizes
</code>

Next, unpack the datasets.

<code>
tar xzf 3dgenome_data.tar.gz
</code>

You should now have 13 files:

<ul>
<li> mm10.chrom.sizes</li>
<li>make_plots.sh</li>
<li>fastq/</li>
	<ul>
		<li>dCTCF/</li>
			<ul>
    		<li>SRR14256290_1.fastq</li>
        <li>SRR14256290_2.fastq</li>
			</ul>
    <li>ddCTCF/</li>
      <ul>
				<li>SRR14256299_1.fastq</li>
        <li>SRR14256299_2.fastq</li>
			</ul>
		</ul>
<li>mm10.DPNII.frag.bed</li>
<li>target.bed</li>
<li>config_hicpro.txt</li>
<li>matrix/</li>
	<ul>
    <li>dCTCF_full.6400.matrix</li>
    <li>ddCTCF_full.6400.matrix</li>
    <li>6400_bins.bed</li>
    <li>dCTCF_full.40000.matrix</li>
    <li>40000_bins.bed</li>
	</ul>
</ul>

The script <code>make_plots.sh</code> is a replacement for a buggy one in HiC-Pro.

The fastq folder contains Hi-C sequenced reads organized by sample. The name <code>dCTCF</code> corresponds to a single deleted CTCF site, while the name <code>ddCTCF</code> corresponds to the CTCF double site deletion.

The file <code>mm10.DPNII.frag.bed</code> contains a list of all of the expected genomic fragments produced by cutting the reference genome with the restriction enzyme DpnII, the one used in these experiments. The file <code>target.bed</code> contains the region that was enriched for (the capture part of capture Hi-C). The file <code>config_hicpro.txt</code> contains settings for running your analysis.

The two <code>.bed</code> files in the matrix raw folders contain the boundaries of each bin produced by breaking the genome into bins of the size indicated in the file name. The <code>.matrix</code> files in the iced folders contain 3 columns corresponding to the bin # of the first part of an interaction, the bin # of the second part of an interaction, and a normalized score associated with that interaction.

<header>
<h4> Setting up HiC-Pro </h4>
</header>

The HiC-Pro software is a complete pipeline for analyzing Hi-C data. It performs the following steps:

<ol>
<li>Mapping reads</li>
<li>Assigning reads to restriction fragments</li>
<li>Quality control</li>
<li>Binning reads into uniform-sized genomic bins</li>
<li>Normalize binned signal</li>
</ol>

To run HiC-Pro, you will first need to build it. Start by activating the conda environment you created at the beginning of the assignment:

<code>
conda activate hicpro
</code>

Next, you will need to clone the HiC-Pro repo. Note that you will be cloning a specific version of the repo rather than the most recent one.

<code>
git clone -b v3.1.0 https://github.com/nservant/HiC-Pro.git
</code>

Finally, you will need to build HiC-Pro. To do this, move into the repo you just cloned and use the following commands:

<code>
make configure PREFIX=${PWD}/../
</code>
<br/><br/>
<code>
make install
</code>

 This will install HiC-Pro into your working directory. Now exit that directory. Finally, you will need to copy the script <code>make_plots.sh</code> into the directory <code>HiC-Pro_3.1.0/scripts/</code>

<header>
<h4> Processing the capture Hi-C data </h4>
</header>

Before running HiC-Pro, you will need to edit the configuration file. There are three lines that you will need to replace <code><YOUR_DIRECTORY></code> with the complete path to the directory you are working in. You can always see the full path with the command <code>pwd</code>. Once you have filled those three entries in, save the config file.

In order to run HiC-Pro, you need to give it the folder where the fastq files are, organized by sample (<code>-i</code>), the name of an output directory to store the results in (<code>-o</code>), and a configuration file (<code>-c</code>). The configuration file was obtained from the repo describing the analysis from the paper and has been updated with the correct paths for various files. You will need to run the program by calling it from <code>HiC-Pro_3.1.0/bin</code> which was created in your working directory when you installed HiC-Pro.

Running the analysis will take several minutes. You will see information about each step as it is performed. Once the analysis has finished running, look around in the output directory. You should have five folders. The <code>bowtie_results</code> directory has all of the mapped read data. The <code>hic_results</code> has all of the data once reads have been assigned to restriction fragments. In this directory, there is a directory <code>pic</code> containing several QC plots. Take a look at the <code>HiCContactRanges</code> and <code>HiCFragment</code> plots (the mapping plots don't tell you anything since only mapped data were selected).

<ul>
	<li><b>What percentage of reads are valid interactions (duplicates do not count as valid)?</b></li>
	<li><b>What constitutes the majority of invalid 3C pairs? What does it actually mean (you may need to dig into the <a href="https://github.com/nservant/HiC-Pro/blob/v3.1.0/doc/MANUAL.md">HiC-Pro manual</a>)?</b></li>
</ul>

</details>
***
