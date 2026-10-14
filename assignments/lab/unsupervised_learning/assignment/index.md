# Unsupervised learning

## Assignment Overview

Today, you will be examining gene expression data from a recent paper by Sahasrabudhe *et al.* entitled ["mTOR signaling contributes to system-driven rhythmic gene expression in mouse liver"](https://pubmed.ncbi.nlm.nih.gov/42490445/). The data are RNA-seq samples taken every 4 hours from mice fed on a nighttime-restricted feeding schedule. Some mice were given rapamycin to inhibit mTOR, a nutrient-sensing kinase responsible for oscillating gene expression patterns due to feeding.

In today's assignment, you will be using unsupervised machine learning methods to explore the data for both quality and bioligical insights, looking at the parsing circadian rhythm genes from genes driven by feeding behaviors.

<br>


## Data

The data you will be utilizing today are transformed read counts (log10 transformed with a pseudocount of one) from three conditions (WT_Ctrl: wild-type control, WT_Rapa: wild-type rapamycin-treated). Each condition consists of 6 time points, sampled every 4 hours and each timepoint consists of 3 replicates. You will use both the complete WT_Ctrl dataset including all of the replicates and datasets with the replicates averaged for WT_Ctrl and WT_Rapa. You also have a list of gene names corresponding to dataset rows and sample names corresponding to the replicate dataset columns.

<br>


## Exercises

There are three exercises in this assignment:

1. Use principle component analysis (PCA) to investigate the variability and relationship between samples and replicates
2. Filter and normalize data and then use PCA to look at the relationship of genes across timepoints
3. Use the PCA results to organize and phase the expressions data followed by k-means clustering to find perturbation patterns from rapamycin treatment.

<br>


## Submission

Before you begin, create a `week5` folder in your `QBB2026-answers` repo. You will be expected to turn in 7 files for this assignment.

1. A python script with comments explaining your code
2. Replicate PCA plot from step 1.3
3. Filtered replicate PCA plot from step 2.1
4. Gene PCA plot from step 2.3
5. Ordered gene expression heatmap from step 3.1
6. Phased gene expression heatmap from step 3.2
7. Clustered gene expression heatmap from step 3.4

<br>


### Exercise 1: Replicate PCA analysis

#### **Step 1.1**: Acquiring the data

To begin, you should download the gzipped file containing all of the data for this assignment. You can do this with the following commands:

```
wget https://github.com/bxlab/cmdb-quantbio/raw/refs/heads/main/assignments/lab/unsupervised_learning/extra_data/unsupervised_learning_data.tgz
tar xzf unsupervised_learning_data.tgz
```

You should have 6 files unpacked:

1. WT_NF_Ctrl_rep_counts.tsv
2. WT_NF_Ctrl_counts.tsv
3. WT_NF_Rapa_counts.tsv
4. gene_names.txt
5. samples.txt
6. unsupervised_framework.py

The first three files are the data you will be working with. The two text files contain labels of the data columns and rows. The last file is a python script containing code to load the data and three helper functions for you to use in the assignmente.

#### **Step 1.2**: Write the PCA analysis (2 points)

To perform PCA, you will use the function `PCA` from the `sklearn.decomposition` module to create a PCA model object. This object can then be trained on your data and provide transformed coordinates for your data. You will then need to plot the first two principal components of your PCA. For this section, you will be using the replicate data as you should see how the timepoints and replicates relate to each other.

 - Create a function that takes in rep data, does a PCA transformation, and plots the first two principal components.
 - Check out documentation for the PCA function https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
 - the PCA should use sample replicates as the x-axis (samples) and gene expression as the y-axis (features)
 - you may need to transpose the matrix, i.e. flip the axes, which you can do using `.T` at the end of your data matrix variable
 - in the PCA function, include the argument `random_state=42`. This will make sure you get the same answer each time
 - to train fit the PCA model and get back the principal components in one command, use the PCA method `.fit_transform`
 - plot the first and second PCs as a scatterplot
 - color the points by the sampling time (`c=times, cmap='Set1'`)
 - size the points by the replicate (`s=reps*15`)
 - it may also be useful to have an argument for a plot title so you can label each plot correctly

#### **Step 1.3**: Plot the first two PCs (0.5 points)

Use the function you wrote to plot the first two PCs for the replicate data. How similar do the replicates look to each other? Would you expect them to be more or less similar than samples from different time points?

Save this plot to turn in.

### Exercise 2: Data filtering and more PCA

#### **Step 2.1**: Filtering data by mean expression and scaled standard deviation (1 point)

Because you don't want to look at every gene, let's limit the data to only genes with sufficient expression and high-enough variation across time. To do this, first find the mean and standard deviation of data for each sample using `numpy.mean` and `numpy.std` and the argument `axis=1` to apply these metrics to each row.

You should now have two 1D arrays of data equal in length to the number of samples (rows), one for means and one for standard deviations. To filter the data using these, you can make them into boolean arrays by putting them into a conditional statement. For example, we only want to keep genes with a mean expression about 0.85 and a scaled standard deviation (i.e. std / mean) greater than 0.115.   

You can use the boolean arrays just like an index array for the data matrix to return only values for which the boolean array is True.

```
filter = value > cutoff
filtered_data = data[:, filter]
```

Filter all three datasets. Once you have the filtered data, use your PCA function to plot the first to PCs of the filtered replicates. Do they look any different?

Save this plot to turn in.

#### **Step 2.2**: Write a function to perform PCA for each gene (1 point)

Next, let's see how the genes you've selected relate to each other. To do this, you will use PCA to reduce the number of dimensions from 6 (the number of timepoints) to 2 (and easy number for plotting). Like in step 1.2, write a function that can find and plot the first two PCs of the WT_NF_Ctrl data. In order to give the plot some context, also accept an argument for specifying a color for each gene via a 1D array.

 - use genes as the x-axis (samples) and timepoints as the y-axis (features)
 - to train fit the PCA model and get back the principal components in one command, use the PCA method `.fit_transform`
 - in the PCA function, include the argument `random_state=42`. This will make sure you get the same answer as us
 - return the first two sets of PC coordinates as you will need them later

#### **Step 2.3**: Normalize data and plot the gene PCA (1 point)

Before using the function you wrote in step 2.2 to plot your filtered gene data for WT_NF_Ctrl, you need to normalize the gene expression across time for each gene. To do this, use the provided `normlize_data` function which returns a normalized matrix along the time axis. Since you've restricted genes to those that have some variation over time, let's identify points by when their peak expression is. This information can then be passed to the color argument. Now you're ready to plot.

 - since you will use them later, normalize both the WT_NF_Ctrl and WT_NF_Rapa data
 - to help understand what the data mean, color the points by their time of peak expression, which you can get using the syntax `peak = numpy.argmax(norm_nf_data, axis=1)`

Save this plot to turn in.

### Exercise 3: Transform and cluster data (1 point)

If everything has gone correctly, you should see that your PCA plot looks like a circle with a relatively smooth color transition around the circle, almost like a clock. If we make the assumption that that is what is being represented, then the angle of each gene's point from the origin should be its expression timing.

#### **Step 3.1**: Order expression times and plot data (1 point)

Let's see what happens when you plot the data using the angle to order the genes. To do this, calculate the angle from the first two PC coordinates (remember that you returned the PC coordinate data from your plotting function) and use the provided `plot_hm` function to visualize the ordered expression data.

 - to find the angles you can use the function `numpy.arctan2(y_coords, x_coords)`
 - the `plot_hm` function has an argument `orderby` which will do the sorting based on an array of continuous or discrete values

#### **Step 3.2**: Phase expression times and plot data (1 point)

Now that you have a way of identify the expression timing offset, you can align the data so each gene is phased to the same timing. This will be important for seeing the different kinds of perturbations in the next step. To phase the data, you can use the provided function `align_times` which takes the data and angle arrays. Once data are aligned, plot the aligned WT_NF_Ctrl data to see that everything looks lined up.

 - align both the WT_NF_Ctrl and WT_NF_Rapa data using the angles calculated from the WT_NF_Ctrl data
 - you only need to plot the WT_NF_Ctrl data

Save this plot to turn in.

#### **Step 3.3**: Cluster data (1.5 points)

Now that the data are aligned, you can finally explore how mTor repression affects the gene expression patterns of cyclicly-expressed genes. using the `KMeans` function from the `sklearn.cluster` module, divide the combined aligned WT_NF_Ctrl and aligned WT_NF_Rapa data into 6 groups.

 - to do k-means clustering, use the function KMeans from sklearn.cluster (see documentation https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html)
 - use `random_state=42` as an argument to KMeans to ensure your clustering is the same as ours
 - use the method `.fit_predict` to fit the model and get predicted labels in a single step
 - to combine the data for predicting groups, you can pass the data in the format `numpy.concatenate((aligned_WT_NF_Ctrl, aligned_WT_NF_Rapa), axis=1)` to the `.fit_predict` function

#### **Step 3.4**: Plot clustered data (1 point) 

Plot the data using these groups to order the data and annotate the groups using the `plot_hm` function provided.

 - the `plot_hm` function accepts multiple positional arguments. Pass both the aligned WT_NF_Ctrl and WT_NF_Rapa matrices to be plotted
 - the `plot_hm` function has an argument `orderby` which will do the sorting based on an array of continuous or discrete values
 - the `plot_hm` function also has an argument `showgroups` which, when set to True, will annotate the groups on the plot

Save this plot to turn in.

<br>

## Grading

1. Function for running PCA on replicates **(2 points)**
2. Plotting replicate PCs **(0.5 points)**
3. Filtering data **(1 point)**
4. Function for running PCA on genes **(1 point)**
5. Normalized data plot of gene PCA **(1 point)**
6. Ordered expression times gene heatmap **(1 point)**
7. Phased expression times gene heatmap **(1 point)**
8. Kmeans clustering **(1.5 points)**
9. Clustered gene heatmap **(1 point)**

**Total Points: 10**

<br>

## Advanced Exercise:

To investigate what the different clusters might mean, use the group labels to extract gene names for each cluster. Don't forget that to do this, the gene names need to filtered in the same way as the genes themselves. To explore what the gene sets may represent, perform a gene ontology enrichment analysis using the website (https://geneontology.org/)[https://geneontology.org/]. Remember that these data are from mouse. You can explore enrichment in three different categories, biological process, molecular function, and cellular component. 
