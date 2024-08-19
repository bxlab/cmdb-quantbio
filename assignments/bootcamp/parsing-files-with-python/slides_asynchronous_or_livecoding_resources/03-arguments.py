#!/usr/bin/env python3

# sys module provides access to command line arguments
import sys

# lists store a sequence of values
print( sys.argv )

# len() returns an integer, + concatenates strings
num_args = len(sys.argv)
print( "There are " + str(num_args) + " arguments" )

# for iterates through lists
i = 0
for arg in sys.argv:
    print( "The " + str(i) + "th argument is " + sys.argv[i] )
    i = i + 1
