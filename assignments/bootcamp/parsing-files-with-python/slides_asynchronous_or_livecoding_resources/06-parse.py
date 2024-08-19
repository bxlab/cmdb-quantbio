#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1] )

for my_line in my_file:
    # in can also check for presence
    if "#" in my_line:
        continue
    # lists store a sequence of values
    fields = my_line.split("\t")
    if fields[2] == "gene":
        print( fields )

my_file.close()
