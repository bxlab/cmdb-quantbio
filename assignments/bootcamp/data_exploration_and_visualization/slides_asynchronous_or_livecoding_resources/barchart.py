#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

genes_per_chrom = np.genfromtxt("genes_per_chrom.txt", 
                                dtype = None,
                                encoding = None,
                                names = ["gene_count", "chrom_name"])
                                
chrom_lengths = np.genfromtxt("https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv", dtype = None, names = ["chrom_name", "chrom_length"])[0:25]

fig, ax = plt.subplots() # create a figure and axes

ax.scatter(chrom_lengths['chrom_length'], genes_per_chrom['gene_count'])

for i, txt in enumerate(genes_per_chrom['chrom_name']): 
  ax.annotate(txt, (chrom_lengths['chrom_length'][i], genes_per_chrom['gene_count'][i]))

ax.set_ylabel("Number of genes/transcripts/exons in GTF")
ax.set_xlabel("Length of chrom. (bp)")

plt.show()





