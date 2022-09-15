#!/usr/bin/env python3

import numpy as np
import sys


#========================#
# Set sequences to align #
#========================#

# In the homework assignment, you'll reading sequences from a
# FASTA file. For the live-coding, we'll define them
# explicitly in the script.

sequence1 = 'TGTTACGG'
sequence2 = 'GGTTGACTA'


#===========================================================#
# Read in match, mismatch, and gap scores from command line #
#===========================================================#

# We can use `sys.argv` to read arguments in from the command
# line, stored in a list. The first thing in the list is
# always the name of the script.



# The rest of the elements of the list are any arguments
# passed when we run the script in the command line



# Assuming three arguments: match score, mismatch score, and
# gap penalty, store these arguments as variables we can use
# in our script.




#============================#
# Playing around with arrays #
#============================#

# Numpy is a Python package built around matrices (arrays).
# You can think of an array as a (potentially)
# multi-dimensional list. In fact, a 2D numpy array is
# essentially a list of lists. We can create a numpy
# array using the `np.array()` function.



# If we want to know what the dimensions of our array
# are, we can check the `.shape` attribute



# Just like we can index lists, we can also index numpy
# arrays. When we index a numpy array, we first index
# the rows, then the columns.



# Like with lists, numpy arrays also support item assignment.



# We can also loop through numpy arrays. When we loop through
# a 2D array, we loop through the rows



# Because each row is just a 1D numpy array, we can also loop
# through the row itself



# We can also use the range() function to loop through each
# value in the array




#=====================#
# Initialize F-matrix #
#=====================#

# The first thing we need to do is create an empty F-matrix.
# The number of rows should be equal to the length of
# sequence1 plus one (to allow for leading gaps). Similarly,
# the number of columns should be equal to the length of
# sequence2 plus one.



# Now we need to fill in the values in the first row and
# first column, based on the gap penalty. Let's fill in the
# first column.



# Now fill in the first row




#=======================#
# Populate the F-matrix #
#=======================#

# Now that we've filled in the first row and column, we need
# to go row-by-row, and within each row go column-by-column,
# calculating the scores for the three possible alignments
# and storing the maximum score




#====================#
# Print the F-matrix #
#====================#




#=========================#
# A primer on while loops #
#=========================#

# While loops are a useful tool in Python (and most other languages)
# that allows you to continue doing a process until some condition
# stops being met. If we want, we can have them mimic a for loop:

# Let's make a for loop that prints out the integers from 0 to 10
# (non-inclusive)


# We can do the exact same thing with a while loop


# While loops can also let us easily go backwards

