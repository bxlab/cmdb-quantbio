# Assignment 8: Bulk RNA-Seq
Assignment Date: Friday, Nov. 11, 2022 <br>
Due Date: Friday, Nov. 18, 2022 <br>

## Lecture and Live-Coding

**Slides** are available here: To be linked

**Master Live-coding script** is available here: `cmdb-quantbio/assignments/lab/machine_learning/slides_asynchronous_or_livecoding_resources/parse_xpresso_predictions.py`. You can copy this to your current directory. Alternatively, you can download it to your current directory with the following curl command: 

`curl https://raw.githubusercontent.com/bxlab/cmdb-quantbio/main/assignments/lab/machine_learning/slides_asynchronous_or_livecoding_resources/parse_xpresso_predictions.py --output parse_xpresso_predictions.py`

This file shows you how we built the script, step-by-step. As such, only the last code block actually needs to be run (we've commented out the "development" code for you). That said, it's possible we don't reach the last code block during class, so if you want to run the code as we had it in class, you may need to un-comment one of the earlier code blocks, and comment the last one.

**Processed Data File for Live-coding** is available here: `~/cmdb-quantbio/assignments/lab/machine_learning/slides_asynchronous_or_livecoding_resources/xpresso_predictions_human.txt`. You can copy this to your current directory. Alternatively, you can download it to your current directory with the following curl command:

`curl https://raw.githubusercontent.com/bxlab/cmdb-quantbio/assignments/lab/machine_learning/slides_asynchronous_or_livecoding_resources/xpresso_predictions_human.txt --output xpresso_predictions_human.txt`

This is the fully processed xpresso data that you'll be using as input to your `parse_xpresso_predictions.py` script. Remember, in class we started with the unprocessed `416685-1.xlsx` excel spreadsheet, and processed into a format that we could more easily used in python. 

## Assignment Overview

Today's assignment is **OPTIONAL**. However, if you do choose to do the assignment, you can have your grade for *this* assignment replace your grade for another of the assignments *of your choice*. That said, we do have some recommendations:
1. If you feel like you are up to speed on all of the material in the course, thus far, we recommend that you do this assignment, and we think that you'll get a lot out of it.
2. If you feel like you're a little bit behind, but the content of this assignment seems exciting to you, we recommend that you do the assignment, and replace your grade for one of the previous assignments that you struggled with. 
3. If you feel like you're a little bit behind, but are not interested in the content of the assignment, you're welcome to take this week to catch up on some of your previous assignments.

With that out of the way: today's assignment comes in two parts. In the first exercise, you'll be following a tutorial to build, from (relative) scratch your own Convolutional Neural Network (CNN) to interpret written numerals. In the second exercise, you'll be using a pre-trained expression prediction model to predict gene expression status in the xpresso data you used during the live-coding exercise. You'll then be comparing the results of running this model with the results that you got during the live-coding exercise.

### Part 1: Building your own CNN

#### Data



