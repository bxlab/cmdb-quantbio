# Assignment 9: Machine Learning 
Assignment Date: Friday, Nov. 11, 2022 <br>
Due Date: Friday, Nov. 18, 2022 <br><br>

## Lecture and Live-Coding

**Slides** are available here: To be linked <br><br>

**Master Live-coding script** is available here: `cmdb-quantbio/assignments/lab/machine_learning/slides_asynchronous_or_livecoding_resources/parse_xpresso_predictions.py`. You can copy this to your current directory. Alternatively, you can download it to your current directory with the following `curl` command: 

`curl https://raw.githubusercontent.com/bxlab/cmdb-quantbio/main/assignments/lab/machine_learning/slides_asynchronous_or_livecoding_resources/parse_xpresso_predictions.py --output parse_xpresso_predictions.py` 

This file shows you how we built the script, step-by-step. As such, only the last code block actually needs to be run (we've commented out the "development" code for you). That said, it's possible we don't reach the last code block during class, so if you want to run the code as we had it in class, you may need to un-comment one of the earlier code blocks, and comment the last one. <br><br>

**Processed Data File for Live-coding** is available here: `~/cmdb-quantbio/assignments/lab/machine_learning/slides_asynchronous_or_livecoding_resources/xpresso_predictions_human.txt`. You can copy this to your current directory. Alternatively, you can download it to your current directory with the following `curl` command:

`curl https://raw.githubusercontent.com/bxlab/cmdb-quantbio/assignments/lab/machine_learning/slides_asynchronous_or_livecoding_resources/xpresso_predictions_human.txt --output xpresso_predictions_human.txt`

This is the fully processed xpresso data that you'll be using as input to your `parse_xpresso_predictions.py` script. Remember, in class we started with the unprocessed `416685-1.xlsx` excel spreadsheet, and processed into a format that we could more easily used in python. <br><br>

## Assignment Overview

Today's assignment is **OPTIONAL**. However, if you do choose to do the assignment, you can have your grade for *this* assignment replace your grade for another of the assignments *of your choice*. That said, we do have some recommendations:

1. If you feel like you are up to speed on all of the material in the course, thus far, we recommend that you do this assignment, and we think that you'll get a lot out of it.
2. If you feel like you're a little bit behind, but the content of this assignment seems exciting to you, we recommend that you do the assignment, and replace your grade for one of the previous assignments that you struggled with. 
3. If you feel like you're a little bit behind, but are not interested in the content of the assignment, you're welcome to take this week to catch up on some of your previous assignments.

With that out of the way: today's assignment comes in two parts. In the first exercise, you'll be following a tutorial to build, from (relative) scratch your own Convolutional Neural Network (CNN) to interpret pictures of different kinds of clothing. In the second exercise, you'll be using a pre-trained expression prediction model to predict gene expression status in the xpresso data you used during the live-coding exercise. You'll then be comparing the results of running this model with the results that you got during the live-coding exercise. <br><br>

## Data

For this assignment you'll need the Python package `tensorflow` that allows you to build, train, and use all sorts of neural networks. We've put together a conda environment that will allow you to use `tensorflow` and several related libraries. To create this conda environment on your computer, you'll need this `.yml` file: `~/cmdb-quantbio/assignments/lab/machine_learning/extra_data/tf_env.yml`. You can copy this to your current directory. Alternatively you can download it to your current directory with the following `curl` command:

`curl https://raw.githubusercontent.com/bxlab/cmdb-quantbio/assignments/lab/machine_learning/extra_data/tf_env.yml --output tf_env.yml`

You can then use this `.yml` file to create a conda environment using the following command:

`conda env create -f tf_env.yml`

This will create a conda environment called `tf_env`, which you can now activate. <br><br>

Additionally, for the second exercise, you'll need some additional from the [Xpresso site](https://xpresso.gs.washington.edu/data/) we looked at during class. Download and unpack the `Xpresso-predict.zip` file. From the unpacked directory, you'll need the following files: `xpresso_predict.py`, `input_fasta/humman_promoters.fa.gz` and `pretrained_models/K562_trainepoch.11-0.4917.h5*`. <br><br>

## Assignment

### Part 1: Building your own CNN

Follow [this tutorial](https://www.tensorflow.org/tutorials/keras/classification) to build and train a CNN that predicts the kind of clothing shown in an input image.

For this exercise, you'll submit the script where you build and train your model, as well as some plots that show your model's permformance (up to you how you want to demonstrate this). Also, in your README, feel free to include any questions you had about the exercise or about NNs in general, and we're happy to address them in our feedback! <br><br>

### Part 2: Using the Xpresso expression prediction model

Using the `xpresso_predict.py` script and the trained model that you downloaded, predict gene expression from the fasta file of human promoters. You can use the `--help` argument with `xpresso_predict.py` to determine what command you need to run. <br><br>

After running the prediction model, compare these results to the results you got from class. Generate a plot that shows the correlation between the predicted expression values you got during the live-coding exercise to the expression values you got from running the Xpresso model. **NOTE**: you may need to ensure that the genes are in the same order between the predicted output here, and predicted expression from the live-coding.

Compare the predictions you got from this model to what was reported in the paper (the txt file we exported from excel and worked with in python in the interactive lecture). In your README, answer the following questions:

1. Are the predictions different, or the same? Are there any trends that you notice looking at the data?
2. Looking at the `xpresso_predict.py` script you downloaded, why do you think the predictions may differ? <br><br>

## Submission

For this assignment, submit the following:

1. A README file including any questions *you* had about this exercise, as well as the command you ran for Part 2 and your answers for the questions in Part 2
2. Your script for Part 1
3. At least one plot demonstrating your model's performance from Part 1
4. Your plotting script for Part 2
5. The correlation plot for Part 2