#!/usr/bin/env python

import sys

import numpy

def main():
    bg = numpy.loadtxt(sys.argv[1], dtype=numpy.dtype([
        ('chr', 'S30'), ('start', int), ('end', int), ('score', float)]))
    bg = bg[numpy.where(bg['chr'] == b"chr17")]
    total = numpy.sum((bg['end'] - bg['start']) * bg['score'])
    total /= numpy.sum(bg['end'] - bg['start'])
    bg['score'] /= total
    output = open(sys.argv[2], mode='w')
    for i in range(bg.shape[0]):
        print(f"{bg['chr'][i].decode('utf8')}\t{bg['start'][i]}\t{bg['end'][i]}\t{bg['score'][i]}",
              file=output)
    output.close()

main()