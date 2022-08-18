#!/usr/bin/env python

"""
Generates violin plots of the Iris dataset.

Also covers:
- Titles
- Axis labels
- xticks
"""

#============================================================#
# IMPORTS

import matplotlib.pyplot as plt
plt.style.use('ggplot')
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

fig, ax = plt.subplots()                 # Open a blank canvas
ax.set_title("Sepal size by species")    # Add a title to the top of the figure
ax.violinplot(                           # Create a violin plot
	[x[Y==0],x[Y==1],x[Y==2]],           # ...of the values for each species
	positions=range(len(species))        # ...at the following positions.
	)
ax.set_xticks(range(len(species)))       # Specify which tick marks to sho.. 
ax.set_xticklabels(species)              # ..and what the labels at the ticks should say
ax.set_xlabel("Species")                 # Specify the label for the x-axis
ax.set_ylabel(labels[0])                 # Specify the label for the y-axis as the first label
fig.savefig("clean_violinplot.png")      # Save the image
plt.close(fig)                           # Close the canvas
