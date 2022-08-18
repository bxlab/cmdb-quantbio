#!/usr/bin/env python

"""
Generates pie chart of the species in the Iris dataset.

Also covers:
- figsize
- bincount
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

fig, ax = plt.subplots(figsize=[7,7])        # Create a blank canvas, 7inches wide by 7inches tall.
                                             ##Keeping it square looks best for a pie chart
ax.set_title("Proportion of each species")   # Add title to the top of the plot
pie = ax.pie(np.bincount(Y), labels=species) #'bincount' counts how many times each integer shows up. Y has values from [0,1,2] 
fig.savefig("clean_pie.png")                 # Save the figure 
plt.close(fig)                               # Close the canvas