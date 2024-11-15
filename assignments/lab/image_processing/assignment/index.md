# Image Processing

## Assignment Overview

The goal of today's lab is to learn how to load images in Python, manipulate values, segment, and score regions of interest across multiple channels. In order to do this, you will be making use of a series if fluorescent microscopy images produced as part of [this study](https://www.nature.com/articles/s41597-021-00944-5) looking at the effects of siRNA-based knock-down of a large array of genes on both RNA production and DNA replication in HeLa cells. The aim of this assignment is to look at the effects of several gene knock-downs compared to a control at single-cell resolution in an automated fashion.

Code from the live-coding session can be downloaded [here](https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/assignments/lab/image_processing/slides_asynchronous_or_livecoding_resources/live_coding.py)

**Important**
Please commit and push as you finish each part. This will really help us gauge the where people are being challenged and better help you all and future students.

<br>

## Data

You will be provided with 24 single-channel images of cells (one color in each image). There are 4 different samples (knocked-down genes), each with 2 separate fields (positions in the sample well). Each sample-field has 3 channel images, one for  DAPI (a stain for DNA), one for Alexa568 (antibody for PCNA, a protein associated with DNA replication), and one for EU-Alexa488 (labeled uracil, marking nascent RNA).

To get the data, download them from [here](https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/assignments/lab/image_processing/extra_data/image_data.tar.gz). You can use the following commands to download and unpack the data:

```bash
wget hhttps://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/assignments/lab/image_processing/extra_data/image_data.tar.gz
tar -xzf image_data.tar.gz
```

<br>

## Exercises

There are three exercises in this assignment:

1. Load the image data into numpy arrays (~20 minutes)
2. Segment and filter each image to label individual nuclei (~1 hour)
3. Find mean signal for each nuclei and plot the results  (~1 hour)

> **_Important:_** Note that in any code snippets shown in the assignment, variable names are examples only and you will need to use the variable names from your own code when applying the functions or operations!

<br>

## Submission

Before you begin, create a `week10` folder in your `QBB2024-answers` repo. You will be expected to turn in 4 files for this assignment.

1. A Python script for loading images, segmenting nuclei, and outputting nuclei signals
2. An R script for plotting the nuclei signals as violin plots
3. A plot with 3 panels of violin plots or 3 separate violin plots
4. A README.md file containing answers to questions from sections 3.1 and 3.2

<br>

### Exercise 1: Loading the image data

Using a `for` loop, load the channel images into python with the `imageio.v3.imread` function, like in the live-coding. You will need to loop through the list of gene names, using the same code to load each gene's data. You may also use a `for` loop for going through each channel and even one for each field if you want.

For each gene/field, you should put all of the channels into the same array such that your image data will end up having a shape of (X, Y, 3) (width, height, # channels). The data type of the image data is numpy.uint16. This will give you a total of 8 3-colored image arrays.

<br>

### Exercise 2: Identifying individual cells

In this section, you will be using the segmentation function that was used in the live-coding to identify cells, then filtering them to remove outliers.

> **_Important:_** If you do the image segmenation and signal measurements all together inside a loop, it will make it easier than having to keep track of each mask, label array, etc.

#### **Step 2.1** For each image, create a binary mask from the DAPI channel

In order to segment the iamges, you will need to create a mask image indicating nucleus vs. background. If you open one of the DAPI images in `Preview`, you will not see anything but a dark field. In order to actually see the data, you need to adjust the color, sliding the upper limit far to the left. This should tell you that the images are rather dim. However, you aren't doing this by eye, you will be using the data to guide you.

To create the mask for each image, use the mean value of an image's DAPI channel as a cutoff, giving you a binary array with `True` for values greater than or equal to the cutoff and `False` for values less than the cutoff.

#### **Step 2.2** Find labels for each image based on the DAPI mask from step 2.1
  
Using the `find_labels` function from the live-coding, create a label map from your DAPI mask. You should get an array of shape (X, Y) with zeros in all non-cell positions and positive numbers denoting pixels belonging to a given nucleus.

You can use `matplotlib.pyplot.imshow(label_array); matplot.pyplot.show()` to see if your labeling makes sense (assuming your label array has the name `label_array`). If you do this, you may want to lower the value of the background with a command like `label_array[numpy.where(label_array == 0)] -= 50`. Be sure to make a copy of your array before you do this, as you will need the background to be zero in a future step.

#### **Step 2.3** Filter out labeled outliers based on size

Because there is always noise in microscope imaging, you may get small dots that are labeled as nuclei. In addition, overlapping cells may result in multiple nuclei labeled as belonging to a single cell. To remove these, you will use size to identify and eliminate them. You can use the `filter_by_size` function from the live coding for this step.

First, filter out any labeled object that is smaller than 100 pixels, as these should be random spots and junk in the image, not nuclei. Remember to set the upper limit high enough not to remove anything large.

Second, using the already filtered label array, find the sizes of each item. To do this, you can use the numpy function `bincount`. This will return essentially a histogram, an array of counts for the number of times each interger appears in the array. When you pass you label array to this function, you will need to flatten it, as `bincount` only accepts 1D arrays. You can easily do this with the array method `ravel()`, e.g. passing the argument `label_array.ravel()` to `bincount`. You will want to ignore the count of zeros, since these are the background (remember how to select only certain ranges using indexing?).

Finally, filter the labels one more time, this time using the mean size plus or minus the standard deviation of sizes as the lower and upper boundaries. This will give you your final labels

<br>

### Exercise 3: Score the PCNA and nascent RNA signal in each nucleus and plot them

Now that you have a set of labeled individual nuclei, you will be finding the mean signal in each nucleus for the two other labeled features, PCNA and nascent RNA. Using these values, you will be able to plot their distributions and test for differences between conditions.

#### **Step 3.1** Find the mean signal for each nucleus from the PCNA and nascent RNA channels

For each gene knock-down condition, you will be constructing a list of values for each nucleus: mean PCNA signal, mean nascent RNA signal, and the log2-transformed ratio of mean nascent RNA signal to mean PCNA signal. Combine data from both fields for a given gene. In the end, you should have three lists of values for each cell, one entry per nucleus in each list.

To find the mean signal for a nucleus, you need to identify which pixel positions correspond to that nucleus. This is where the label array comes in. For example, to find the pixels for nucleus #7, you would use the command `where = numpy.where(label_array == 7)`. The number of labeled nuclei for a given label array will be equal to `numpy.amax(label_array) + 1`. Don't forget to skip label zero, as this corresponds to the background. You will need a `for` loop to go through each label in a given label array.

Write these data to a text file. I suggest using the column format `Gene,nascentRNA,PCA,ratio` as this will make it easier to plot with R.

**What do each of these values mean, based on the descriptions of what is being labeled?**

#### **Step 3.2** Plot each set of data in a separate violin plot

For each set of data, nascent RNA signal, PCNA signal, and the log2 ratio, create a violin plot with appropriate axis labels and title.

**Look at the knocked-down gene with the highest ratio and the one with the lowest ratio. Look up what their functions are and explain if you think this result makes sense.**

<br><br>

## Grading

1. Python script **(6 pts total)**:
  * Used a `for` loop to cycle through genes **(0.5 pts)**
  * Combined channels into a single array for each image **(0.5 pts)**
  * Created a binary mask from the DAPI channel **(0.5 pts)**
  * Created mask and label array from DAPI channel **(0.5 pts)**
  * Filtered out small objects **(0.5 pts)**
  * Found object sizes, excluding background **(0.5 pts)**
  * Found correct minimum and maximum size cutoffs for second filtering **(0.5 pts)**
  * Filtered out too small and large objects **(0.5 pts)**
  * Correctly identified pixels for a given label **(0.5 pts)**
  * Found signal for correct 2 channels and log-transformed ratio for each nucleus **(1 pt)**
  * Wrote values to text file **(0.5 pts)**
2. R script **(0.5 pts total)**:
  * Load and plot data in violin plots **(0.5 pts)**
3. Violin plots **(0.5 pts total)**:
  * Created 3 properly labeled violin plots **(0.5 pts)**
4. README.md **(1 pt total)**:
  * Answered question 3.1 **(0.5 pts)**
  * Answerd question 3.2 **(0.5 pts)**

**5. Make 2 commits on your project work by 11/22 (2 pts)**

**Total Points: 10**

<br><br>

## Advanced Excercises (Optional)

### Exercise 4: Test the significance of the differences in conditions

In this experiment, PIM2 served as a control as it was expected to have no effect on either DNA replication or RNA production. Using a two-tailed Student's T-test, determine which, if any, conditions were significantly different from the PIM2 condition for each of the measurments. 

### Exercise 5: Plot the data as color images with all channels shown together

`Matplotlib.pyplot.imshow` has the ability to show color images. The requirement is that the array passed is 3D, with shape width x height x RGB and that the array is of type `numpy.uint8`, i.e. all values are integers from 0-255. The other challenge is that in order to best visualize all three channels together, they need to be independently scaled before displaying the image. This means scaling by the min and max values and then ensuring that values are rounded to the nearest integer and within the 0-255 value constraints.
