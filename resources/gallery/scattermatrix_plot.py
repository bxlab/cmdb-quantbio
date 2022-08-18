#!/usr/bin/env python

"""
Generates scatter matrix plots of the Iris dataset.

Also covers:
- Titles
- Axis labels
- Subplots
- Adjust subplots
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

n = X.shape[1]                         # Store the number of features in X

fig, axes = plt.subplots(              # Create a blank canvas
	nrows=n, ncols=n,                  # Generate enough axes for each variable combination
	figsize=(10,10)                    # Set overall figure size to 10 inches x 10 inches
	) 

for i in range(n):                     # Iterate through the indexes of the features
	for j in range(n):                 # ... and again.
		ax = axes.flat[i*n+j]          # ... ... Select the one in the (i+1)-th row in the (j+1)-th position in that row 
		if i == j:                     # ... ... If the indices match
			ax.hist(X[:,i],            # ... ... ... Create a histogram of the values
			 edgecolor='lightgray')    #             with light gray edges
		else:                          # ... ... Else:
			ax.scatter(X[:,j],X[:,i],  # ... ... ... Create a scatter plot of the j-th vs i-th features. 
			marker='.', color='b')     #             with blue points as markers
		if ax.is_last_row():           # ... ... If the plot is in the bottom row
			ax.set_xlabel(labels[j])   # ... ... ... Add a label to the x-axis
		if ax.is_first_col():          # ... ... If the plot is in the first column
			ax.set_ylabel(labels[i])   # ... ... ... Add a label to the y-axis

fig.subplots_adjust(                   # Adjust the spacing of the subplots, to help make everything fit
    left = 0.075,                      # ... the left edge of the left-most plot will be this percent of the way across the width of the plot
    bottom = 0.05,                     # ... the bottom edge of the bottom-most plot will be this percent of the way up the canvas
    right = 0.98,                      # ... the right edge of the right-most plot will be this percent of the way across the width
    top = 0.98,                        # ... the top edge of the top-most plot will be this percent of the way from the bottom
    wspace = 0.2,                      # ... the space between the plots side-to-side
    hspace = 0.2,                      # ... the space between the plots top-to-bottom
)

fig.savefig("scattermatrix.png")       # Save the figure
plt.close(fig)                            # Close the canvas