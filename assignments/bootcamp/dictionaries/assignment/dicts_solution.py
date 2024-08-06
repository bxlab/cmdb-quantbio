#!/usr/bin/env python

import sys

import numpy


def main():
    attribute_fname, gene_tissue_fname, data_fname, out_fname = sys.argv[1:5]

    # Load gene-tissue pair data
    GT_dict = {}
    for line in open(gene_tissue_fname):
        gene_ID, gene_name, tissue = line.rstrip().split("\t")
        GT_dict[gene_ID] = tissue

    # Load sample IDs for tissue association
    sample_dict = {}
    fs = open(attribute_fname)
    _ = fs.readline()
    for line in fs:
        fields = line.rstrip().split("\t")
        sample = fields[0]
        tissue = fields[6]
        sample_dict.setdefault(tissue, [])
        sample_dict[tissue].append(sample)

    # Load expression data
    genes = []
    data = []
    fs = open(data_fname)
    
    # Skip header lines
    _ = fs.readline()
    _ = fs.readline()
    
    # Save sample IDs for colum  identification
    sample_IDs = fs.readline().rstrip().split("\t")[2:]

    # Convert sample IDs to column indices in sample dict
    for tissue in sample_dict.keys():
        indices = []
        for sample_ID in sample_dict[tissue]:
            # Check that sample_ID is in data sample_IDs
            if sample_ID in sample_IDs:
                indices.append(sample_IDs.index(sample_ID))
        sample_dict[tissue] = numpy.array(indices)
        print(tissue, len(indices))

    # Load data for each gene
    for line in fs:
        fields = line.rstrip().split('\t')
        gene = fields[0]
        # If not a marker gene, skip
        if gene not in GT_dict:
            continue
        # Record gene to keep track of order
        genes.append(gene)
        # Convert expressions to a 1D array
        expr = numpy.array(fields[2:], dtype=float)
        # Get tissue associated with marker gene
        tissue = GT_dict[gene]
        # Get tissue sample indices
        sample_indices = sample_dict[tissue]
        # Keep tissue expressions
        data.append(expr[sample_indices])
    fs.close()

    # Save expressions by tissue
    output = open(out_fname, 'w')
    output.write(f"GeneID\tTissue\tExpr\n")
    for i in range(len(genes)):
        gene = genes[i]
        tissue = GT_dict[gene]
        for j in range(len(data[i])):
            output.write(f"{gene}\t{tissue}\t{data[i][j]}\n")
    output.close()

main()