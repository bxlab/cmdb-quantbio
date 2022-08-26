
import sys

fname = sys.argv[1]
fs = open(fname, 'r')
bed = []
for i, line in enumerate(fs):
    if line.startswith("#"):
        continue
    fields = line.rstrip().split("\t")
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    bed.append([chrom, start, end])
fs.close()
