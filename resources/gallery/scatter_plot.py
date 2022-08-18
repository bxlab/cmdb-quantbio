#!/usr/bin/env python

"""
Generates scatter plots of the Iris dataset.

Also covers:
- colors
- marker styles
- Boolean filtering
- Legend
- Titles
- Axis labels
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

# Extract two features from the data to plot
x = X[:,2]   
y = X[:,3]


fig1, ax1 = plt.subplots()         # Create a blank canvas
ax1.scatter(x,y)                   # Plot x vs y as points
fig1.savefig("first_scatter.png") # Save the figure
plt.close(fig1)                   # Close the canvas


colors = ['darkblue','orange','cyan'] # List some colors to plot the points with
ms = ['o','d','s']                    # Choose some marker styles for plotting, a circle, diamond, and a square.

fig2, ax2 = plt.subplots()              # Create a new blank canvas
ax2.set_title("Petal Size\nby species") # Add a title to the top, spanning two lines
for i in range(len(species)):           # Iterate through indices corresponding to species names
	ax2.scatter(                        # ...Create a scatter plot
		x[Y==i],                        # ... ... of x
		y[Y==i],                        # ... ... vs y
		c=colors[i],                    # ... ... set the color of these points
		marker=ms[i],                   # ... ... and pick the marker style
		label=species[i]                # ... ... finally, label the points of this set with the species name
		)
l = ax2.legend(loc='upper left')        # Add a legend to the top left
ax2.set_xlabel(labels[2])               # label the x-axis
ax2.set_ylabel(labels[3])               # label the y-axis
fig2.savefig("clean_scatter.png")       # Save the figure
plt.close(fig2)                         # Close the canvas
