#!/usr/bin/env python3

import sys

def parse_bed(fname):
    fs = open(fname, 'r')
    bed = []
    field_types = [str, int, int, str, float, str]
    for i, line in enumerate(fs):
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        fieldN = len(fields)
        if fieldN < 3:
            print(f"Line {i} appears malformed", file=sys.stderr)
            continue
        for j in range(min(len(field_types), fieldN)):
            fields[j] = field_types[j](fields[j])
        bed.append(fields)
    fs.close()
    return bed

if __name__ == "__main__":
    fname = sys.argv[1]
    bed = parse_bed(fname)
