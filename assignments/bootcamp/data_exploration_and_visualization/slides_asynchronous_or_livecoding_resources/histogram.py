#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

# np.random.seed(42)
# norm_data = np.random.normal(0, 1, int(1e6))

snp_af = np.genfromtxt("chr21_af.txt")
#print(snp_af.shape)

plt.hist(snp_af, bins = 50)

plt.show()