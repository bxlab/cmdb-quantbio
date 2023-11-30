# QuantLab Week 8 - Data Visualization
Assignment Date: Friday, Dec. 1, 2023 <br>
Due Date: Friday, Dec. 8, 2023 <br>

## Lecture

**Slides** are available here: [Lecture slides](https://www.dropbox.com/scl/fi/tytab80ncq1ia85remlsh/20231201_qblab_dataviz.pptx?rlkey=x08ydaut17vn3x8dha0zdfebo&dl=0)


## Assignment Overview

In today's assignment we will practice the principles of data visualization that we covered in the lecture. Visualizations should accurately and concisely convey a message in an easily understandable way without misleading the audience. Avoid extraneous visual elements.

We will start with some more directed visualizations of the GTEx RNA-seq dataset you analyzed last week, followed by a more open-ended project.

## Data

Re-load the data you downloaded in the previous lab session. These data comprise RNA-seq samples from whole blood from the GTEx Consortium (755 total individuals). They were downloaded directly from the [GTEx portal](https://gtexportal.org/home/downloads/adult-gtex#bulk_tissue_expression) and slightly reformatted to save you some time on tedious data wrangling. If needed, download again from the Dropbox links below:

[Subject-level metadata](https://www.dropbox.com/scl/fi/zidlbn4rlvyv43k022mmn/gtex_metadata.txt?rlkey=j6aidakljr0739tnnzvpbg0gn&dl=0) </br>
[Gene expression matrix](https://www.dropbox.com/scl/fi/7iengpyrevd356dfq53pg/gtex_whole_blood_counts_formatted.txt?rlkey=l5h12cyher33kkzlrwi4qwf8g&dl=0)

Load and normalize the data using the code below:

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

# merge with metadata
full_design_df = pd.concat([counts_df_normed, metadata], axis=1)
```

### Step 1: Visualize characteristics of the RNA-seq data (barplots, line plots, box/violin plots, and histograms)

Generate the following figures. For all figures, label the axes appropriately, provide legends only when necessary, do not place a title or any other elements on the plot itself. 

1. For subject GTEX-113JC, plot the distribution of expression across all 54,592 genes. Use a [log scale](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xscale.html#matplotlib.axes.Axes.set_xscale) if necessary.

2. Plot the expression distribution of the gene MXD4 in males versus females.

3. Plot the number of subjects per age category that fall into each category of the Hardy scale.

4. Plot the median expression of the gene LPXN for each age category, stratified by sex.

### Step 2: Independent data visualization

Working with a partner, select any dataset from the TidyTuesday repository on GitHub: https://github.com/rfordatascience/tidytuesday. 

1. Using Python (`pandas`, `numpy`, `matplotlib`, etc.) these data with your partner, searching for any interesting features or patterns. Jot these down as notes (no need to submit).

2. Choose three aspects of these data that are best represented by three different types of figures (e.g., histogram, bar plot, line plot, heatmap, etc.).

3. Generate these figures with `matplotlib`. As always, label the axes appropriately and avoid extraneous visual elements. Using `matplotlib.pyplot.title`, write a title on your plot that concisely states the message that your figure is attempting to convey.


## Submission

For this assignment you should submit the following:

1. `plotting_step1.py` script to load and analyze data and produce plots from Step 1

2. `plotting_step2.py` script to load and analyze data and produce plots from Step 2

3. The pretty plots produced by the scripts above.

**DO NOT** push any raw data! Only the things we asked for!<br><br>


## Additional Resources

* [Fundamentals of Data Visualization Textbook](https://clauswilke.com/dataviz/)

## Advanced Exercises

We'd like to see if we can identify any broad patterns present in our gene expression data. To explore this, we're going to cluster the data, both by sample as well as by gene.

As a first step, we will need to log2-transform our data:

```
counts_df_normed = np.log2(counts_df_normed + 1)
```

To perform clustering, you'll be using the `dendrogram`, `linkage` and `leaves_list` functions from `scipy`. The documentation for SciPy isn't very helpful, but with some quick Googling you can find examples of how to use both of these tools.

1. Using `linkage` and `leaves_list`, cluster the filtered and log2 transformed gene expression data matrix for both genes and samples based on their patterns of expression (so both the rows and columns of the matrix). You will find the numpy transpose functionality useful in order to cluster the columns.

2. Plot a heatmap of the clustered gene expression data.

3. Using the dendrogram function, create a dendrogram relating the samples to one another.
