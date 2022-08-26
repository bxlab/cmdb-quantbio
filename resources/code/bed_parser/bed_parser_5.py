#!/usr/bin/env python3

import sys

"""Usage: python bedParser.py <BED_FILE>"""

def parse_bed(fname):
    # Make sure the file exists
    try:
        fs = open(fname, 'r')
    except:
        raise FileNotFoundError("That file doesn’t appear to exist")
    bed = []
    malformed = 0
    # Define data types for all 12 possible fields
    field_types = [str, int, int, str, float, str,
                   int, int, str, int, str, str]
    for i, line in enumerate(fs):
        # Check for comment lines
        if line.startswith("#"):
            continue
        # Split line into fields, removing newline characters
        fields = line.rstrip().split("\t")
        # Check that the number of fields is appropriate
        fieldN = len(fields)
        if fieldN != 12 and (fieldN < 3 or fieldN > 10):
            malformed += 1
            continue
        # Safely try converting fields to correct data type
        try:
            for j in range(fieldN):
                fields[j] = field_types[j](fields[j])
        except:
            malformed += 1
            continue
        # Convert itemRGB if needed and check range
        if fieldN > 8:
            if fields[8] == "0":
                fields[8] = 0
            else:
                # Make sure itemRGB has exactly 3 values
                fields[8] = fields[8].split(",")
                if len(fields[8]) != 3:
                    malformed += 1
                    continue
                # Ensure itemRGB is all ints
                try:
                    for j in range(3):
                        fields[8][j] = int(fields[8][j])
                except:
                    malformed += 1
                    continue
        # If needed, convert and check blockSizes and blockStarts
        if fieldN > 10:
            fields[10] = fields[10].rstrip(',').split(',')
            fields[11] = fields[11].rstrip(',').split(',')
            # Make sure the lists aren't too short and can be ints
            try:
                for j in range(fields[9]):
                    fields[10][j] = int(fields[10][j])
                    fields[11][j] = int(fields[11][j])
            except:
                malformed += 1
                continue
            # Check block lengths
            if len(fields[10]) != fields[9] or len(fields[11]) != fields[9]:
                malformed += 1
                continue
        bed.append(fields)
    # If there were any malformed entries, report them
    if malformed > 0:
        print(f"There were {malformed} malformed entries", file=sys.stderr)
    fs.close()
    return bed

if __name__ == "__main__":
    fname = sys.argv[1]
    bed = parse_bed(fname)
