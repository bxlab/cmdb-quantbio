#!/usr/bin/env python

import sys

import numpy


def main():
    # Load GTEx filename and destination for results
    gene_fname, min_max_fname, gene_tissue_fname = sys.argv[1:4]

    # Open file
    fs = open(gene_fname, 'r')
    
    # Skip first two lines
    _ = fs.readline()
    _ = fs.readline()

    # Read tissue names
    tissues = fs.readline().rstrip().split("\t")[2:]

    # Create lists to hold results in
    gene_IDs = []
    gene_names = []
    data = []

    # Read in genes
    for line in fs:
        fields = line.rstrip().split("\t")
        gene_IDs.append(fields[0])
        gene_names.append(fields[1])
        data.append(fields[2:])

    # close file
    fs.close()

    # convert data into arrays
    gene_IDs = numpy.array(gene_IDs)
    gene_names = numpy.array(gene_names)
    data = numpy.array(data, float)

    print(numpy.median(data), numpy.mean(data))

    # Log-transform expression
    data = numpy.log2(1 + data)

    print(numpy.median(data), numpy.mean(data))

    # Find expression stats for each gene
    gene_mins = numpy.amin(data, axis=1)
    gene_maxes = numpy.amax(data, axis=1)

    # Save first 5000 mins, maxes
    output = open(min_max_fname, 'w')
    output.write("Gene\tMin\tMax\n")
    for i in range(5000):
        output.write(f"{gene_names[i]}\t{gene_mins[i]}\t{gene_maxes[i]}\n")
    output.close()

    # Find difference between highest and second highest expression
    sorted = numpy.sort(data, axis=1)
    diffs = sorted[:, -1] - sorted[:, -2]

    # Identify highly tissue-specific genes and their tissues
    high_tissue = numpy.zeros_like(data)
    high_tissue[numpy.arange(data.shape[0]), numpy.argmax(data, axis=1)] = 1
    marker_genes = (diffs > 10).reshape(-1, 1) * high_tissue# (data > 10)
    high_GTs = numpy.where(marker_genes == 1)

    # Print gene-tissue pairs
    output = open(gene_tissue_fname, 'w')
    for i in range(high_GTs[0].shape[0]):
        output.write(f"{gene_IDs[high_GTs[0][i]]}\t{gene_names[high_GTs[0][i]]}\t{tissues[high_GTs[1][i]]}\n")
    output.close()

main()