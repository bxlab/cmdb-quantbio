# Data Visualization

## Assignment Overview

In today's assignment you will practice the principles of data visualization that were covered in the lecture. Visualizations should accurately and concisely convey a message in an easily understandable way without misleading the audience. Avoid extraneous visual elements.

You will start with some more directed visualizations of the GTEx RNA-seq dataset you analyzed last week, followed by a more open-ended project.

## Data

Re-load the data you downloaded in the previous lab session. These data comprise RNA-seq samples from whole blood from the GTEx Consortium (755 total individuals). As a reminder, they were downloaded directly from the [GTEx portal](https://gtexportal.org/home/downloads/adult-gtex#bulk_tissue_expression) and slightly reformatted to save you some time on tedious data wrangling. If needed, download again from the Dropbox links below:

[Subject-level metadata](https://www.dropbox.com/scl/fi/zidlbn4rlvyv43k022mmn/gtex_metadata.txt?rlkey=j6aidakljr0739tnnzvpbg0gn&dl=0)<br>
[Gene expression matrix](https://www.dropbox.com/scl/fi/7iengpyrevd356dfq53pg/gtex_whole_blood_counts_formatted.txt?rlkey=l5h12cyher33kkzlrwi4qwf8g&dl=0)<br><br>

## Exercises

### Excercise 1: Visualize characteristics of the RNA-seq data

In the first exercise, you'll be exploring some aspects of the GTEx whole blood data and generating plots that communicate the observed patterns in an easy-to-interpret manner. You will be producing four figures. Each should be saved as it's own separate file and uploaded as part of the assignment.

For all figures, label the axes appropriately, provide legends only when necessary, do not place a title or any other elements on the plot itself. **You will be graded on proper labelling.**

Create a `plotting_exercise1.py` script now for this exercise.<br><br>

#### **Step 1.0**: Load the data

In your `plotting_exercise1.py` script, load and normalize the GTEx data using the code below:

```
import numpy as np
import pandas as pd
from pydeseq2 import preprocessing
from matplotlib import pyplot as plt

# read in data
counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

# read in metadata
metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)

# normalize
counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]

# log
counts_df_logged = np.log2(counts_df_normed + 1)

# merge with metadata
full_design_df = pd.concat([counts_df_logged, metadata], axis=1)
```
<br>

#### **Step 1.1**: Distribution of expression across genes

For subject GTEX-113JC, plot the distribution of expression (logged normalized counts) across all genes, excluding any genes with 0 counts. Upload this figure for the assignment.<br><br>

#### **Step 1.2**: Expression of a single gene between sexes

For the gene MXD4, plot the distribution of expression (logged normalized counts) in males versus females. Upload this figure for the assignment.<br><br>

#### **Step 1.3**: Distribution of subject ages

Plot the number of subjects in each age category. Upload this figure for the assignment.<br><br>

#### **Step 1.4**: Sex-stratified expression with age

For the gene LPXN, plot the median expression (logged normalized counts) over time (i.e. in each age category), stratified by sex. Upload this figure for the assignment.<br><br>

### Exercise 2: Independent data visualization

We highly encourage you to do this exercise with a partner.

With your partner, select any data set from the TidyTuesday repository on GitHub: [https://github.com/rfordatascience/tidytuesday](https://github.com/rfordatascience/tidytuesday). In this exercise, you'll be exploring interesting patterns in the data set that you choose.

Create a `plotting_exercise2.py` script now for this exercise.<br><br>

#### **Step 2.1**: Initial exploration

Using Python (`pandas`, `numpy`, `matplotlib`, etc.) explore these data with your partner, searching for any interesting features or patterns. Jot down any interesting patterns you observe as notes (no need to submit). For each feature/pattern you observe, think about what kind of plot would best communicate that feature/pattern.<br><br>

#### **Step 2.2**: Visualization

Choose three aspects/patterns of these data that are best represented by three different types of plots (e.g., histogram, bar plot, line plot, heatmap, etc.).

Generate these figures using `matplotlib`. As always, label the axes appropriately and avoid extraneous visual elements. For each plot, add a title that concisely states the message that your figure is attempting to convey.<br><br>

## Submission

For this assignment you should submit the following:

1. `plotting_exercise1.py` script to load and analyze data and produce plots from Exercise 1 **(3 points total)**
    * Code to load in and normalize data **(1 point)**
    * Code to create plot for Step 1.1 **(0.5 point)**
    * Code to create plot for Step 1.2 **(0.5 point)**
    * Code to create plot for Step 1.3 **(0.5 point)**
    * Code to create plot for Step 1.4 **(0.5 point)**
2. `plotting_exercise2.py` script to load and analyze data and produce plots from Exercise 2 **(3.5 points total)**
    * Code to load in data **(1 point)**
    * Code for initial data exploration **(1 point)**
    * Code to produce first plot **(0.5 point)**
    * Code to produce second plot **(0.5 point)**
    * Code to produce third plot **(0.5 point)**
3. Pretty plots produced by the scripts above. **(3.5 points total)**
    * Plot for Step 1.1 **(0.5 point)**
    * Plot for Step 1.2 **(0.5 point)**
    * Plot for Step 1.3 **(0.5 point)**
    * Plot for Step 1.4 **(0.5 point)**
    * First plot for Exercise 2 **(0.5 point)**
    * Second plot for Exercise 2 **(0.5 point)**
    * Third plot for Exercise 2 **(0.5 point)**

**Total Points: 10**

**DO NOT** push any raw data! Only the things we asked for!<br><br>

## Advanced Exercises

### Advanced Exercise 1

For each age range, plot the proportion of samples within each group of the Hardy scale. This should be a single panel figure.<br><br>

### Advanced Exercise 2

We'd like to see if we can identify any broad patterns present in our gene expression data. To explore this, we're going to cluster the data, both by sample as well as by gene.

To perform clustering, you'll be using the `dendrogram`, `linkage` and `leaves_list` functions from `scipy`. The documentation for SciPy isn't very helpful, but with some quick Googling you can find examples of how to use both of these tools.

1. Using `linkage` and `leaves_list`, cluster the filtered and log2 transformed gene expression data matrix for both genes and samples based on their patterns of expression (so both the rows and columns of the matrix). You will find the numpy transpose functionality useful in order to cluster the columns.

2. Plot a heatmap of the clustered gene expression data.

3. Using the dendrogram function, create a dendrogram relating the samples to one another.<br><br>

## Additional Resources

* [Fundamentals of Data Visualization Textbook](https://clauswilke.com/dataviz/)<br><br>
