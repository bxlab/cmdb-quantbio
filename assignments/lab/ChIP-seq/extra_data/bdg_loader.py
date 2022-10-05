#!/usr/bin/env python

import sys

import numpy

def load_data(fname):
    # Load in the bedgraph with columns chr, start, end, and score
    bg = numpy.loadtxt(fname, usecols=(1,2,3), dtype=numpy.dtype([
        ('start', int), ('end', int), ('score', float)]))
    # Find the lowest coordinate
    start = bg['start'][0]
    # Create an array at bp resolution
    array = numpy.zeros(bg['end'][-1] - start, dtype=numpy.dtype(
        [('pos', int), ('score', float)]))
    # For each line in the bedgraph, set positions across range with score
    for i in range(bg.shape[0]):
        array['score'][bg['start'][i] - start:bg['end'][i] - start] = bg['score'][i]
    # Set coordinates in array
    array['pos'] = numpy.arange(start, start + array.shape[0])
    # Bin array into 100bp bins
    hist = numpy.histogram(array['pos'], weights=array['score'], bins=(array.shape[0] // 100))
    # Create final data dict
    X = (hist[1][1:] + hist[1][:-1]) / 2
    Y = hist[0]
    data = {'X': X, 'Y':Y}
    return data