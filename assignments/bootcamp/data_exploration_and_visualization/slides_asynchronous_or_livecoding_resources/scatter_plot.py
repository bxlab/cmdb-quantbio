#!/usr/bin/env python

import matplotlib.pyplot as plt

x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]

x2 = [2, 4, 6]
y2 = [8, 64, 216]

# create a figure and axes
fig, ax = plt.subplots(nrows = 2)
#(len(ax))
ax[0].plot(x, y, label = "x^2")
ax[1].plot(x2, y2, label = "x^3")
ax[0].legend()
plt.savefig("lineplot.png")
plt.close(fig)
