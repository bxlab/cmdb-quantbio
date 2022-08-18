#!/usr/bin/env python

"""
Generates plots of the Iris dataset.

Also covers:
- np.argsort
- Titles
- Axis labels
- Line styles
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

# Grab the first two columns of X
x = X[:,0]
y = X[:,1]


fig1, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 6)) # Open a blank canvas, 8 inches x 6 inches
                                                         # with two subplots, next to each other
# Plot the columns with default parameters in one subplot
ax1.plot(x,y)

# Sort x and y based on x, and plot these against
# each other
indexes = np.argsort(x)  # 'argsort' returns the indices of the 
                         # values from smallest to largest
x = x[indexes]           # We can use these indices to sort
y = y[indexes]           # x and y.

ax2.plot(x,y)

# Save the figure as a .png file, and close the canvas
fig1.savefig("first_plot.png")
plt.close(fig1)


# Plot a nice clean plot, with labels
fig2, ax3 = plt.subplots(figsize=(8,6))        # Open canvas
ax3.set_title("Sepal dimensions")              # Will appear above the plot
plot = ax3.plot(x, y, linestyle=':',color='k') # Plots x vs y with a dotted black line
ax3.set_xlabel(labels[0])                      # Label the x axis
ax3.set_ylabel(labels[1])                      # Label the y axis
fig2.savefig("clean_plot.png")                 # Save
plt.close(fig2)                                # Close canvas