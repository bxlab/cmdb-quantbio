# Image Classification with PyTorch

## Overview

Computer vision involves tasks such as image classification, object detection, instance segmentation, and action classification.
The Image Processing lab showed how edges (features) can be detected using kernels (filters).
With enough good data, [deep learning](https://pytorch.org) can provide a way to discover more features and filters.
In this assignment, you will investigate how much data is needed to accurately classify [hand-written digits](https://en.wikipedia.org/wiki/MNIST_database) and how well [models](https://github.com/bioimage-io/use-cases/tree/main/case3-devtools) can segment images from the [Human Protein Atlas](https://www.proteinatlas.org).

## Exercises

### 1. PyTorch -- How much data is needed?

- 1a. What is the accuracy using 60,000 images? 30,000? 6,000? 3,000? 600?
- 1b. How do the weights looks different when trained with 60,000 vs 600?

### 2. BioImage -- How well does segmentation work?

- 2a. Explain the parameters provided to watershed(). What happens when you remove mask? What happens when you remove markers?
- 2b. Compare how well segmentation works for cyto, endo, mito, and nucl. What types of samples work well? Not as well?

## Submit

Submit two Jupyter Notebooks with your code, output, and answers.
Start by making a copy of the analysis.ipynb and add your work at the bottom of the notebook.
Format your answers in Markdown blocks and include enough code and output to show your work.

PyTorch

- Question 1a answer (1 pt)
- Question 1b answer (1 pt)
- pytorch.ipynb well organized (3 pt)

BioImage

- Question 2a answer (1 pt)
- Question 2b answer (1 pt)
- bioimage.ipynb well organized (3 pt)
