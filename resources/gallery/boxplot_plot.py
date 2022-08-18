#!/usr/bin/env python

"""
Generates boxplots of the Iris dataset.

Also covers:
- Titles
- Axis labels
"""

#============================================================#
# IMPORTS

import matplotlib.pyplot as plt
from sklearn import datasets
import numpy as np

#============================================================#
# LOAD SOME DATA TO USE

iris = datasets.load_iris()
X = iris.data
Y = iris.target
labels = ["Sepal Length", "Sepal Width", "Petal Length", "Petal Width",]
species = ["Setosa", "Versicolour", "Virginica"]

#============================================================#
# CREATE SOME PLOTS

# Extract the first feature from the dataset
x = X[:,0]

fig, ax = plt.subplots()              # Open a blank canvas
ax.set_title("Sepal size by species") # Add a title to the top
ax.boxplot(                           # Create a violin plot
	[x[Y==0],x[Y==1],x[Y==2]],        # ...of the values for each species
	labels=species,                   # ...with the species names as labels
	)
ax.set_xlabel("Species")              # Label the x-axis
ax.set_ylabel(labels[0])              # Label the y-axis with the first feature name
fig.savefig("clean_boxplot.png")      # Save the plot
plt.close(fig)                        # Close the canvas
