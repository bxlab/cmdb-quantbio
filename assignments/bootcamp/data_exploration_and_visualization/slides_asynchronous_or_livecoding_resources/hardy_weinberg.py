#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from vcfParser import parse_vcf

af_list = []
gt00_list = []
gt01_list = []
gt11_list = []

fname = sys.argv[1]
vcf = parse_vcf(fname)

for i in range(1, len(vcf)):
    snp = vcf[i]
    af_list.append(snp[7]["AF"])
    gt00_count = 0
    gt01_count = 0
    gt11_count = 0
    for j in range(9, len(snp)):
        if (snp[j] == "0|0"):
            gt00_count += 1
        elif (snp[j] == "0|1" or snp[j] == "1|0"):
            gt01_count += 1
        elif (snp[j] == "1|1"):
            gt11_count += 1
    gt00_list.append(gt00_count / 2548)
    gt01_list.append(gt01_count / 2548)
    gt11_list.append(gt11_count / 2548)

# create a figure and axes
fig, ax = plt.subplots()
ax.scatter(af_list, gt00_list, label = "Hom. Ref.")
ax.scatter(af_list, gt01_list, label = "Het.")
ax.scatter(af_list, gt11_list, label = "Hom. Alt.")
ax.set_ylabel("Genotype Freq.")
ax.set_xlabel("Allele Freq.")
ax.legend()

# add theoretical expectations
x = np.arange(0.0, 1.0, 0.01)
hom_ref = x**2
het = 2 * x * (1 - x)
hom_alt = (1 - x)**2
ax.plot(x, hom_ref, lw = 2, color = "black", linestyle = "dashed")
ax.plot(x, het, lw = 2, color = "black", linestyle = "dashed")
ax.plot(x, hom_alt, lw = 2, color = "black", linestyle = "dashed")

plt.show()

            
    
    
    
    



