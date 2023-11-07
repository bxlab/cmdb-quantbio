#!/usr/bin/env python

import sys

import matplotlib.pyplot as plt


def main():
    # Load file names from command line
    normal_fname, tumor_fname, out_fname = sys.argv[1:4]

    # Load data from files
    normal = load_data(normal_fname)
    tumor = load_data(tumor_fname)

    # Find reads that appear more than once in datasets
    normal_set = set()
    normal_multi = set()
    for i in range(len(normal)):
        if normal[i][0] not in normal_set:
            normal_set.add(normal[i][0])
        else:
            normal_multi.add(normal[i][0])
    normal_single = normal_set.difference(normal_multi)

    tumor_set = set()
    tumor_multi = set()
    for i in range(len(tumor)):
        if tumor[i][0] not in tumor_set:
            tumor_set.add(tumor[i][0])
        else:
            tumor_multi.add(tumor[i][0])
    tumor_single = tumor_set.difference(tumor_multi)

    # Print statistics about unique vs multimapping reads
    print(f"Normal unique reads: {len(normal_single)} ({len(normal_single) / len(normal_set) * 100}) %")
    print(f"Normal multi reads: {len(normal_multi)} ({len(normal_multi) / len(normal_set) * 100}) %")
    print(f"Tumor unique reads: {len(tumor_single)} ({len(tumor_single) / len(tumor_set) * 100}) %")
    print(f"tumor multi reads: {len(tumor_multi)} ({len(tumor_multi) / len(tumor_set) * 100}) %")

    # Parse data into unique and multi-mapping reads for plotting

    # Plot 4 panels, normal and tumor, unique and multi-mapping

# Create a reusable plotting function
#def plot_data(ax, X, Y, label):

# Create a reusable data loading function. Record read name, actual size, and mapped size
def load_data(fname):
    data = []
    for line in open(fname):
        line = line.rstrip().split()
        data.append([
            line[0], int(line[1]), int(line[3]) - int(line[2])])
    return  data

main()





















