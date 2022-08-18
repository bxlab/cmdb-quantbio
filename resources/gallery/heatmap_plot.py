#!/usr/bin/env python

"""
Generates heatmap of the Iris dataset.

Also covers:
- Title
- imshow
- colormap
- grid
- label ticks
- colorbar
- adjust subplots
- normalizing
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

# Start by normalizing the data so each feature can fit on the
##same scale. Basically we are calculating the Z-score for 
##each value based on only the feature it belongs to.
X = (X-np.average(X,axis=0))/np.std(X,axis=0)

# Pull out the value with the greatest magnitude, to set the
##scale
m = np.max(np.abs(X))

# Make the actual plot
fig, ax = plt.subplots(figsize=(8, 6))          # Open a blank canvas, 8 inches x 6 inches
ax.set_title("Heatmap of Iris characteristics") # Add a title to the top
im = ax.pcolor(                              # Treat the values like pixel intensities in a picture
	X,                                       # ... Using X as the values
	cmap="RdBu",                             # ... Use the Red-white-blue colormap to assign colors to your pixel values
	vmin=-1*m,                               # ... Set the lowest value to show on the scale
	vmax=m,                                  # ... Set the highest value to show on the scale. Since we are using a 'diverging' colormap, these should match.
	)

ax.grid(False)                      # Turn of the grid lines (a feature added automatically by ggplot)
ax.set_xticks(                      # Edit the xticks being shown
	np.arange(0.5, X.shape[1]+0.5), # ... use the values centered on each column of pixels
	)
ax.set_xticklabels(                 # Label the ticks
	labels,                         # ... at position which correspond to the indices of our labels
	rotation=50,                    # ... and rotate the labels 50 degrees counter-clockwise
	)
ax.set_yticks([])                   # Edit the ticks on the y-axis to show....NOTHING

cbar = fig.colorbar(im, ax=ax)      # Add a bar to the right side of the plot which shows the scale correlating the colors to the pixel values

fig.subplots_adjust( # Adjust the spacing of the subplots, to help make everything fit
    left = 0.05,     # ... the left edge of the left-most plot will be this percent of the way across the width of the plot
    bottom = 0.15,   # ... the bottom edge of the bottom-most plot will be this percent of the way up the canvas
    right = 1.0,     # ... the right edge of the right-most plot will be this percent of the way across the width
    top = 0.95,      # ... the top edge of the top-most plot will be this percent of the way from the bottom
)

fig.savefig("clean_heatmap.png") # Save the image
plt.close(fig) # Close the canvas
