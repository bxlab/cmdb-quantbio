#!/usr/bin/env python

###############################################################################

# We're downloading the excel data from here: https://xpresso.gs.washington.edu/data/
# We added comment/hashes/pound signs to the comment lines on the excel spreadsheet, and saved to a txt tab delimited file

"""
import sys

for i,line in enumerate(open(sys.argv[1])):
   if not line.startswith("#"):
       fields = line.strip('\r\n').split('\t')
       if i < 10:
           print(fields)
           quit()
"""

###############################################################################

# We tried to read in the txt file, but found that everything had quotes around it, so we weren't filtering out the header lines.... what can we do to remove these quotes? Something in excel? nothing obvious. What does google say? seems microsoft specific involving notepad. Something in python? Let's strip the quotes away!

"""
import sys

for i,line in enumerate(open(sys.argv[1])):
   if not line.strip('"').startswith("#"):
       fields = line.strip('"\r\n').split('\t')
       if i < 10:
           print(fields)
           quit()
"""

###############################################################################

# Now let's save the gene names into a list and the predicted expression from the k562 model into a list and then the true expression from the k562 cells into a list (so 3 lists and we want column 2 (index 1), 5 (index 4), and ... what the .. which column for the true k562 expression?)

# We can look at the second page in the excel spreadsheet, and we see that we want the E123.

# But how do we know which column is E123? Ugh. Let's save the header header line into an array that we can search to find the index rather than counting. Went back to the excel spreadsheet and added a second # to that

"""
import numpy as np
import sys

for i, line in enumerate(open(sys.argv[1])):
   if line.strip('"').startswith("##"):
       header = line.strip('"\r\n').split('\t')
       print(header)
       k562_obs_idx = np.where(header == "E123")[0]
       print(k562_obs_idx)
   elif not line.strip('"').startswith("#"):
       fields = line.strip('"\r\n').split('\t')
       if i < 10:
           print(fields)
           quit()
"""

###############################################################################

# When we try to find the index, it's blank, and that seems to be because ...., any ideas?

"""
import numpy as np
import sys

for i, line in enumerate(open(sys.argv[1])):
   if line.strip('"').startswith("##"):
       header = np.array(line.strip('"\r\n').split('\t'))
       print(header)
       k562_obs_idx = np.where(header == "E123")[0]
       print(k562_obs_idx)
   elif not line.strip('"').startswith("#"):
       fields = line.strip('"\r\n').split('\t')
       if i < 10:
           print(fields)
           quit()
"""

###############################################################################

# We now know that we want column 61 (index 60) for the E123, and we've saved that into a variable for ourselves.
# let's save our 3 lists of data
# Do we want to do anything for the model predictions and observations data when we save it? Like looking at the printed fields, do we need to do data conversions maybe? -- yes

"""
import numpy as np
import sys

gene_names = []
k562_model_predictions = []
k562_observations = []
for line in open(sys.argv[1]):
    if line.strip('"').startswith("##"):
        header = np.array(line.strip('"\r\n').split('\t'))
        k562_obs_idx = np.where(header == "E123")[0]
    elif not line.strip('"').startswith("#"):
        fields = line.strip('"\r\n').split('\t')
        gene_names.append(fields[1])
        k562_model_predictions.append(float(fields[4]))
        k562_observations.append(float(fields[k562_obs_idx]))
print(len(gene_names))
print(len(k562_model_predictions))
print(len(k562_observations))

print(k562_model_predictions[0:5])
print(k562_observations[0:5])
"""

###############################################################################

# Ugh there's a pretty unhelpful error that shows up, but we can handle it -- the error
#Traceback (most recent call last):
#  File "parse_xpresso_predictions.py", line 82, in <module>
#    k562_observations.append(float(fields[k562_obs_idx]))
#TypeError: only integer scalar arrays can be converted to a scalar index

# after googling "indexing a list python TypeError: only integer scalar arrays can be converted to a scalar index"
#https://stackoverflow.com/questions/50997928/typeerror-only-integer-scalar-arrays-can-be-converted-to-a-scalar-index-with-1d

#either turn fields into a numpy array, or do another [0] around the idx variable because that variable was in a list-like-object when we printed it

"""
import numpy as np
import sys

gene_names = []
k562_model_predictions = []
k562_observations = []
for line in open(sys.argv[1]):
   if line.strip('"').startswith("##"):
       header = np.array(line.strip('"\r\n').split('\t'))
       k562_obs_idx = np.where(header == "E123")[0]
   elif not line.strip('"').startswith("#"):
       fields = line.strip('"\r\n').split('\t')
       gene_names.append(fields[1])
       k562_model_predictions.append(float(fields[4]))
       k562_observations.append(float(fields[k562_obs_idx[0]]))
print(len(gene_names))
print(len(k562_model_predictions))
print(len(k562_observations))

print(k562_model_predictions[0:5])
print(k562_observations[0:5])
"""

###############################################################################

# Now we want to make a scatter plot comparing the observations to the predictions
# predictions on x
# true on y

"""
import numpy as np
import sys
import matplotlib.pyplot as plt

gene_names = []
k562_model_predictions = []
k562_observations = []
for line in open(sys.argv[1]):
    if line.strip('"').startswith("##"):
        header = np.array(line.strip('"\r\n').split('\t'))
        k562_obs_idx = np.where(header == "E123")[0]
    elif not line.strip('"').startswith("#"):
        fields = line.strip('"\r\n').split('\t')
        gene_names.append(fields[1])
        k562_model_predictions.append(float(fields[4]))
        k562_observations.append(float(fields[k562_obs_idx[0]]))

fig, ax = plt.subplots()
ax.scatter(k562_model_predictions, k562_observations, alpha=1, color="blue", s=0.25)
ax.set_xlabel("Predicted K562 expression level,\n10-fold cross-validated")
ax.set_ylabel("K562 expression level (log10)")
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
plt.tight_layout()
plt.show()
"""

###############################################################################

# Add the trendline and r^2, but dude, they did abline(0,1, red), not a trendline, and they used the cor function to find the correlation/r^2
# Our manual for our r^2: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html
# Our manual for adding text to the plot: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html
# our manual for adding a trendline: https://www.python-graph-gallery.com/scatterplot-with-regression-fit-in-matplotlib

"""
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

gene_names = []
k562_model_predictions = []
k562_observations = []
for line in open(sys.argv[1]):
   if line.strip('"').startswith("##"):
       header = np.array(line.strip('"\r\n').split('\t'))
       k562_obs_idx = np.where(header == "E123")[0]
   elif not line.strip('"').startswith("#"):
       fields = line.strip('"\r\n').split('\t')
       gene_names.append(fields[1])
       k562_model_predictions.append(float(fields[4]))
       k562_observations.append(float(fields[k562_obs_idx[0]]))

fig, ax = plt.subplots()
ax.scatter(k562_model_predictions, k562_observations, alpha=1, color="blue", s=0.25)
ax.set_xlabel("Predicted K562 expression level,\n10-fold cross-validated")
ax.set_ylabel("K562 expression level (log10)")
corr = pearsonr(k562_model_predictions, k562_observations)
ax.text(.5, 3.75, "r^2 = " + str(round(corr.statistic**2,2)) + "\nn = " + str(len(k562_observations)))
line_xs = np.linspace(min(min(k562_model_predictions), min(k562_observations)), max(max(k562_model_predictions), max(k562_observations)), 100)
line_ys = 0 + 1 * line_xs
ax.plot(line_xs, line_ys, color = "maroon")
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
plt.tight_layout()
plt.show()"""

###############################################################################

# If we compare this to the figure in the paper, what do we notice is different on the axes and the "trendline"?

"""
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

gene_names = []
k562_model_predictions = []
k562_observations = []
for line in open(sys.argv[1]):
   if line.strip('"').startswith("##"):
       header = np.array(line.strip('"\r\n').split('\t'))
       k562_obs_idx = np.where(header == "E123")[0]
   elif not line.strip('"').startswith("#"):
       fields = line.strip('"\r\n').split('\t')
       gene_names.append(fields[1])
       k562_model_predictions.append(float(fields[4]))
       k562_observations.append(float(fields[k562_obs_idx[0]]))

fig, ax = plt.subplots()
ax.scatter(k562_model_predictions, k562_observations, alpha=1, color="blue", s=0.25)
ax.set_xlabel("Predicted K562 expression level,\n10-fold cross-validated")
ax.set_ylabel("K562 expression level (log10)")
corr = pearsonr(k562_model_predictions, k562_observations)
ax.text(.5, 3.75, "r^2 = " + str(round(corr.statistic**2,2)) + "\nn = " + str(len(k562_observations)))
line_xs = np.linspace(max(min(k562_model_predictions), min(k562_observations)), min(max(k562_model_predictions), max(k562_observations)), 100)
line_ys = 0 + 1 * line_xs
ax.plot(line_xs, line_ys, color = "maroon")
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
plt.tight_layout()
plt.show()
"""

###############################################################################

#"finally" (but not really if we have time or someone points out the real final difference...), we need to make some of the genes red, and labeled
# source for gene list: https://github.com/vagarwal87/Xpresso/blob/master/Fig3_S3/Fig3ABCDEF_S3ABC.R line 214 though. 

"""
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

gene_names = []
k562_model_predictions = []
k562_observations = []
descriptions = []
for line in open(sys.argv[1]):
    if line.strip('"').startswith("##"):
        header = np.array(line.strip('"\r\n').split('\t'))
        k562_obs_idx = np.where(header == "E123")[0]
    elif not line.strip('"').startswith("#"):
        fields = line.strip('"\r\n').split('\t')
        gene_names.append(fields[1])
        k562_model_predictions.append(float(fields[4]))
        k562_observations.append(float(fields[k562_obs_idx[0]]))
        descriptions.append(fields[2])

genesoi = ["PIM1", "SMYD3", "FADS1", "PRKAR2B", "GATA1", "MYC"]
genesoilocs =[]
for geneoi in genesoi:
    genesoilocs.append(np.where(np.array(gene_names) == geneoi)[0][0])
for i in range(len(descriptions)):
    if "hemoglobin subunit" in descriptions[i]:
        genesoi.append(gene_names[i])
        genesoilocs.append(i)

fig, ax = plt.subplots()
ax.scatter(k562_model_predictions, k562_observations, alpha=1, color="blue", s=0.25)
ax.set_xlabel("Predicted K562 expression level,\n10-fold cross-validated")
ax.set_ylabel("K562 expression level (log10)")
corr = pearsonr(k562_model_predictions, k562_observations)
ax.text(.5, 3.75, "r^2 = " + str(round(corr.statistic**2,2)) + "\nn = " + str(len(k562_observations)))
line_xs = np.linspace(max(min(k562_model_predictions), min(k562_observations)), min(max(k562_model_predictions), max(k562_observations)), 100)
line_ys = 0 + 1 * line_xs
ax.plot(line_xs, line_ys, color = "maroon")
for i, geneoi in enumerate(genesoi):
    ax.text(k562_model_predictions[genesoilocs[i]], k562_observations[genesoilocs[i]], geneoi, color="maroon", fontweight="demi")
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
plt.tight_layout()
plt.show()
"""

###############################################################################

# OK the final final difference if we have time or someone points it out
# so what about the cloud like look that the R graph from the paper has??? It's because they use an R function called smoothScatter and we can approximate something close with Python, but let me tell you, it's annoying and we're going to follow a matplotlib walkthrough
# googled "smoothscatter ggplot equivalent for matplotlib"
# found this: https://www.inwt-statistics.com/read-blog/smoothscatter-with-ggplot2-513.html
# knew I could use a stat density 2d equivalent for matplotlib
# googled this: "stat density 2d matplotlib" and found this: https://www.python-graph-gallery.com/2d-density-plot/
# scrolling until I found the walkthrough I wanted
# matplotlib walkthrough: https://www.python-graph-gallery.com/85-density-plot-with-matplotlib
# picking our colormap: https://matplotlib.org/stable/tutorials/colors/colormaps.html

###############################################################################

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, kde #note that the walkthrough says kde not gaussian_kde, but I couldn't figure out how to change it successfully quickly

gene_names = []
k562_model_predictions = []
k562_observations = []
descriptions = []
for line in open(sys.argv[1]):
    if line.strip('"').startswith("##"):
        header = np.array(line.strip('"\r\n').split('\t'))
        k562_obs_idx = np.where(header == "E123")[0]
    elif not line.strip('"').startswith("#"):
        fields = line.strip('"\r\n').split('\t')
        gene_names.append(fields[1])
        k562_model_predictions.append(float(fields[4]))
        k562_observations.append(float(fields[k562_obs_idx[0]]))
        descriptions.append(fields[2])

genesoi = ["PIM1", "SMYD3", "FADS1", "PRKAR2B", "GATA1", "MYC"]
genesoilocs =[]
for geneoi in genesoi:
    genesoilocs.append(np.where(np.array(gene_names) == geneoi)[0][0])
for i in range(len(descriptions)):
    if "hemoglobin subunit" in descriptions[i]:
        genesoi.append(gene_names[i])
        genesoilocs.append(i)

fig, ax = plt.subplots()
x, y = np.array(k562_model_predictions), np.array(k562_observations)
nbins = 300
k = kde.gaussian_kde([x, y])
xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto', cmap=plt.cm.Blues)
ax.set_xlabel("Predicted K562 expression level,\n10-fold cross-validated")
ax.set_ylabel("K562 expression level (log10)")
corr = pearsonr(k562_model_predictions, k562_observations)
ax.text(.5, 3.75, "r^2 = " + str(round(corr.statistic**2,2)) + "\nn = " + str(len(k562_observations)))
line_xs = np.linspace(max(min(k562_model_predictions), min(k562_observations)), min(max(k562_model_predictions), max(k562_observations)), 100)
line_ys = 0 + 1 * line_xs
ax.plot(line_xs, line_ys, color = "maroon")
for i, geneoi in enumerate(genesoi):
    ax.text(k562_model_predictions[genesoilocs[i]], k562_observations[genesoilocs[i]], geneoi, color="maroon", fontweight="demi")
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
plt.tight_layout()
plt.show()
