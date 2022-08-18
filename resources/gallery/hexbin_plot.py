#!/usr/bin/env python

"""
Generates scatterplot and density plot of the Iris dataset.

Also covers:
- Title vs SupTitle
- normalizing
- subplot
- hexbin
- colormap
- axis limits
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

# Start by normalizing the data
X = (X-np.average(X,axis=0))/np.std(X,axis=0)

# Pull out the features we want to plot
x = X[:,0]
y = X[:,2]

# Plot the figure
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 6)) # Open a blank canvas, 8 inches x 6 inches
                                                        # with two subplots, next to each other
fig.suptitle("Distribution of %s vs %s" % (labels[0], labels[2])) # A 'suptitle' is a title for the WHOLE figure, as opposed to just a single plot

ax1.set_title("Scatterplot")              # Add a title to this subplot
ax1.scatter(x, y, color='b', marker='.')  # Create a scatter plot in this canvas

ax2.set_title("Density")                  # Add a title to this subplot
hb = ax2.hexbin(                          # Create a 'hexbin' a plot that uses hexagons to show the density of values
	x,y,                                  # These are the points we will be using
	gridsize=15,                          # The number of hexagons wide to make the image. Fewer bins means more points in each bin
	cmap="YlOrRd"                         # Set the colormap. Here we use the sequential "Yellow-Orange-Red" map.
	) 
#cbar = fig.colorbar(hb, ax=ax2)          # Add a colorbar to the subplot

xmin,xmax = ax2.get_xlim()                # Ask the plot to tell us what it set the axis limits at
ymin,ymax = ax2.get_ylim()                # Ask the plot to tell us what it set the axis limits at
ax1.set_xlim(xmin,xmax)                   # Set the limits of this axis to match that from the hexbin plot
ax1.set_ylim(ymin,ymax)                   # Set the limits of this axis to match that from the hexbin plot



fig.savefig("clean_density.png")   # Save the figure
plt.close(fig)                     # Close the canvas