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

print('Script name:', sys.argv[0])

# The rest of the elements of the list are any arguments
# passed when we run the script in the command line

print('Arguments passed:', sys.argv[1:])

# Assuming three arguments: match score, mismatch score, and
# gap penalty, store these arguments as variables we can use
# in our script.

match_score = float(sys.argv[1])
mismatch_score = float(sys.argv[2])
gap_penalty = float(sys.argv[3]) # Should be negative


#============================#
# Playing around with arrays #
#============================#

# Numpy is a Python package built around matrices (arrays).
# You can think of an array as a (potentially)
# multi-dimensional list. In fact, a 2D numpy array is
# essentially a list of lists. We can create a numpy
# array using the `np.array()` function.

list_of_lists = [[1,2,3],
				 [4,5,6],
				 [7,8,9]]

my_2D_array = np.array(list_of_lists)

print(my_2D_array)

# If we want to know what the dimensions of our array
# are, we can check the `.shape` attribute

print(my_2D_array.shape) # a tuple of (nrows, ncols)

# Just like we can index lists, we can also index numpy
# arrays. When we index a numpy array, we first index
# the rows, then the columns.

my_list = [3,1,4,1,5]

print(my_list[2]) # print element at index 2 in `my_list`
print(my_list[0:3]) # print first three elements of `my_list`

print(my_2D_array[1,2]) # Print `6`
print(my_2D_array[1,:]) # Print the second (i.e. index 1) row
print(my_2D_array[:,1]) # Print the second (i.e. index 1) column
print(my_2D_array[:2,:2]) # Print the upper left 2x2 block

# Like with lists, numpy arrays also support item assignment.

my_list[4] = 6 # Change the `5` to a `6`
print(my_list)

my_2D_array[1,1] = 314 # Change the `5` to a `314`
print(my_2D_array) 

# We can also loop through numpy arrays. When we loop through
# a 2D array, we loop through the rows

for row in my_2D_array:
	print(row)

# Because each row is just a 1D numpy array, we can also loop
# through the row itself

for row in my_2D_array:
	for val in row:
		print(val)

# We can also use the range() function to loop through each
# value in the array

for row_index in range(my_2D_array.shape[0]): # The number of rows is the first value in .shape
	for col_index in range(my_2D_array.shape[1]): # The number of columns is the second value in .shape
		print(my_2D_array[row_index, col_index])


#=====================#
# Initialize F-matrix #
#=====================#

# The first thing we need to do is create an empty F-matrix.
# The number of rows should be equal to the length of
# sequence1 plus one (to allow for leading gaps). Similarly,
# the number of columns should be equal to the length of
# sequence2 plus one.

F_matrix = np.zeros((len(sequence1)+1, len(sequence2)+1))

# Now we need to fill in the values in the first row and
# first column, based on the gap penalty. Let's fill in the
# first column.

for i in range(len(sequence1)+1):
	F_matrix[i,0] = i*gap_penalty

# Now fill in the first row

for j in range(len(sequence2)+1):
	F_matrix[0,j] = j*gap_penalty


#=======================#
# Populate the F-matrix #
#=======================#

# Now that we've filled in the first row and column, we need
# to go row-by-row, and within each row go column-by-column,
# calculating the scores for the three possible alignments
# and storing the maximum score

for i in range(1, len(sequence1)+1): # loop through rows
	for j in range(1, len(sequence2)+1): # loop through columns
		if sequence1[i-1] == sequence2[j-1]: # if sequence1 and sequence2 match at positions i and j, respectively...
			d = F_matrix[i-1, j-1] + match_score
		else: # if sequence1 and sequence2 don't match at those positions...
			d = F_matrix[i-1, j-1] + mismatch_score
		h = F_matrix[i,j-1] + gap_penalty
		v = F_matrix[i-1,j] + gap_penalty

		F_matrix[i,j] = max(d,h,v)


#====================#
# Print the F-matrix #
#====================#

print(F_matrix)


#=========================#
# A primer on while loops #
#=========================#

# While loops are a useful tool in Python (and most other languages)
# that allows you to continue doing a process until some condition
# stops being met. If we want, we can have them mimic a for loop:

# Let's make a for loop that prints out the integers from 0 to 10
# (non-inclusive)
for i in range(0, 10, 1):
	print(i)

# We can do the exact same thing with a while loop
i = 0
while i < 10:
	print(i)
	i += 1 # increment i

# While loops can also let us easily go backwards
i = 9
while i >= 0:
	print(i)
	i -= 1 # decrement i
