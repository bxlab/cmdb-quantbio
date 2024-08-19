#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1] )

# for iterates through files
for my_line in my_file:
    # objects have methods
    my_line = my_line.rstrip("\n")
    print( my_line )

my_file.close()
