
import sys

def parse_bed(fname):
    fs = open(fname, 'r')
    bed = []
    field_types = [str, int, int, str, float, str]
    for i, line in enumerate(fs):
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        for j in range(min(len(field_types), len(fields))):
            fields[j] = field_types[j](fields[j])
        bed.append(fields)
    fs.close()
    return bed

fname = sys.argv[1]
bed = parse_bed(fname)
