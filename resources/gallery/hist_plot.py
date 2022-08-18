#!/usr/bin/env python

"""
Generates histogram of the Iris dataset.

Also covers:
- Boolean filtering
- Titles
- Axis labels
- Alpha transparency
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

fig, ax = plt.subplots(figsize=(8, 6))  # Open a blank canvas, 8 inches x 6 inches
hist = ax.hist(x)                       # Generate a histogram of the data, with defaul settings
fig.savefig("first_hist.png")           # Save the figure
plt.close(fig)                          # Close the canvas

minimum = np.min(x)           # Grab the minimum...
maximum = np.max(x)           # ...and maximum value to set the range

fig2, ax2 = plt.subplots(figsize=(8, 6))          # Open a blank canvas, 8 inches x 6 inches
ax2.set_title("Sepal size by species") # Add a title to the top
for i in range(len(species)):      # Iterate through indices corresponding to the species labels
	n, bins, patches = ax2.hist(   # ... plot a histogram of
		x[Y==i],                   # ... ... only the values corresponding to that species
		bins=30,                   # ... ... Use thirty bars
		range=[minimum,maximum],   # ... ... ranging from the minimum to the maximum
		normed=True,               # ... ... Normalize the bars to frequencies instead of counts
		label=species[i],          # ... ... Assign the corresponding species name as the label
		alpha=0.5,                 # ... ... Make the bars only 50% opaque
		edgecolor='lightgray',     # ... ... Outline the bars with gray
		)                          
l = ax2.legend(loc="upper right")  # Add a legend of with the species labels, to the upper right of the plot.
ax2.set_ylabel("Frequency")        # Label the y-axis
ax2.set_xlabel(labels[0])          # Label the x-axis with the name of the feature we used
fig2.savefig("clean_hist.png")     # Save figure as .png
plt.close(fig2)                    # Close the canvas
