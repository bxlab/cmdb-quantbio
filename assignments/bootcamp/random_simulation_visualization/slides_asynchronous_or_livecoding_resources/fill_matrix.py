#!/usr/bin/env python

import numpy

probabilities = numpy.around(numpy.arange(0, 0.5, 0.1), decimals = 2)
tosses = numpy.arange(10, 40, 10)

new_twodim_arr = numpy.zeros((len(probabilities), len(tosses)))
for i, prob in enumerate(probabilities):
    for j,toss in enumerate(tosses):
        new_twodim_arr[i,j] = prob + toss
print(new_twodim_arr)
