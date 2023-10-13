# QuantLab Week 5 - Variant Calling

Assignment Date: Friday, Oct. 13, 2023

Due Date: Friday, Oct. 27, 2023

## Lecture -- Rajiv McCoy

[Lecture Slides](https://www.dropbox.com/scl/fi/q356m8nya6qv0cd8y5o4y/20231012_qblab_variant_calling.pptx?rlkey=9pltuif66aasxnvmlt0koclo2&dl=0)

## Setting up your environment

There are a bunch of command line tools you'll need for both the live-coding and assignment. You'll be creating a `conda` environment for this week that contains well-behaving versions of all of these tools. Create and activate a new `variants` conda environment for this week using the code below:

```
mamba create -n variants
mamba activate variants
conda config --env --set subdir osx-64
mamba install python=3.7 snpeff=5 freebayes bwa vcflib samtools
```

## Live-coding

As always, before you do anything else, create a `week5` directory in your `qbb2023-answers` directory for this assignment, and create a `README.md` file in that `week5` directory. This is where you will be putting and uploading all of your code/plots/etc. for this assignment.

We would recommend creating a `livecoding` directory within your `week5` directory, and putting your live-coding scripts and data there. You can add the `livecoding` directory to your `.gitignore` so that `git` won't track these files.

Wherever you choose to do your live-coding, you can download the following file:
[https://www.dropbox.com/scl/fi/ye4h2eemcvsomsaccag4h/live_coding_data.tar.gz?rlkey=lgwmjno4ojojxbdwzj4fjamai&dl=0](https://www.dropbox.com/scl/fi/ye4h2eemcvsomsaccag4h/live_coding_data.tar.gz?rlkey=lgwmjno4ojojxbdwzj4fjamai&dl=0)

This file contains the data you'll need for the live-coding exercises. Once downloaded, unzip it with the following command:

```
tar -xvzf live_coding_data.tar.gz
```

<!--
## Homework Assignment

 Complete the homework assignment in your `week4` submission directory in your `qbb2023-answers`.

 [Homework Assignment](../assignments/lab/variant_calling/assignment/index.html)
 -->
