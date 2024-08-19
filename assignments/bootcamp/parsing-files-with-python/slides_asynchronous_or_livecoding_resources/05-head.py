#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1] )

max_lines = 10
if len(sys.argv) > 2:
    max_lines = int(sys.argv[2])

# be careful of off-by-1 errors
i = 0
for my_line in my_file:
    if i >= max_lines:
        break
    my_line = my_line.rstrip("\n")
    print( my_line )
    i = i + 1

my_file.close()
