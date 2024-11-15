#!/usr/bin/env python

import numpy
import scipy
import matplotlib.pyplot as plt
import imageio
import plotly
import plotly.express as px

na = numpy.newaxis

# Load image into numpy array
img = imageio.v3.imread("illum_DAPI.tif")
img = img.astype(numpy.float32) - numpy.amin(img)
img /= numpy.amax(img)

img1 = imageio.v3.imread("illum_RNAcytoNuc.tif")
img1 = img1.astype(numpy.float32) - numpy.amin(img1)
img1 /= numpy.amax(img1)

img2 = imageio.v3.imread("illum_Mito.tif")
img2 = img2.astype(numpy.float32) - numpy.amin(img2)
img2 /= numpy.amax(img2)

rgbimg = numpy.zeros((1080, 1080, 3), numpy.uint16)
for i, name in enumerate(['DAPI', 'RNAcytoNuc', 'Mito']):
    rgbimg[:, :, i] = imageio.v3.imread(f"illum_{name}.tif")

# plt.imshow(rgbimg)
# plt.show()

rgbimg = rgbimg.astype(numpy.float32)
for i in range(3):
    rgbimg[:, :, i] -= numpy.amin(rgbimg[:, :, i])
    rgbimg[:, :, i] /= numpy.amax(rgbimg[:, :, i])
print(numpy.amax(rgbimg.reshape(-1, 3), axis=0))
rgbimg = (numpy.minimum(255, numpy.floor(rgbimg * 256))).astype(numpy.uint8)
plt.imshow(rgbimg)
plt.show()

numpy.mean(img)

mask = img >= 0.035
plt.imshow(mask, cmap='grayscale')
plt.show()

def find_labels(mask):
    # Set initial label
    l = 0
    # Create array to hold labels
    labels = numpy.zeros(mask.shape, numpy.int32)
    # Create list to keep track of label associations
    equivalence = [0]
    # Check upper-left corner
    if mask[0, 0]:
        l += 1
        equivalence.append(l)
        labels[0, 0] = l
    # For each non-zero column in row 0, check back pixel label
    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                # If back pixel has a label, use same label
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[0, y] = l
    # For each non-zero row
    for x in range(1, mask.shape[0]):
        # Check left-most column, up  and up-right pixels
        if mask[x, 0]:
            if mask[x - 1, 0]:
                # If up pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                # If up-right pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[x, 0] = l
        # For each non-zero column except last in nonzero rows, check up, up-right, up-right, up-left, left pixels
        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    # If up pixel has label, use that label
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    # If not up but up-right pixel has label, need to update equivalence table
                    if mask[x - 1, y - 1]:
                        # If up-left pixel has label, relabel up-right equivalence, up-left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x - 1, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x - 1, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    elif mask[x, y - 1]:
                        # If left pixel has label, relabel up-right equivalence, left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    else:
                        # If neither up-left or left pixels are labeled, use up-right equivalence label
                        labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    # If not up, or up-right pixels have labels but up-left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    # If not up, up-right, or up-left pixels have labels but left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    # Otherwise, add new label
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l
        # Check last pixel in row
        if mask[x, -1]:
            if mask[x - 1, -1]:
                # if up pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                # if not up but up-left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                # if not up or up-left but left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x, -2]]
            else:
                # Otherwise, add new label
                l += 1
                equivalence.append(l)
                labels[x, -1] = l
    equivalence = numpy.array(equivalence)
    # Go backwards through all labels
    for i in range(1, len(equivalence))[::-1]:
        # Convert labels to the lowest value in the set associated with a single object
        labels[numpy.where(labels == i)] = equivalence[i]
    # Get set of unique labels
    ulabels = numpy.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[numpy.where(labels == j)] = i
    return labels

labels = find_labels(mask)
sizes = numpy.bincount(labels.ravel())


for i in range(1, numpy.amax(labels)+1):
    where = numpy.where(labels == i)
    if sizes[i] < 100:
        labels[where] = 0


def filter_by_size(labels, minsize, maxsize):
    # Find label sizes
    sizes = numpy.bincount(labels.ravel())
    # Iterate through labels, skipping background
    for i in range(1, sizes.shape[0]):
        # If the number of pixels falls outsize the cutoff range, relabel as background
        if sizes[i] < minsize or sizes[i] > maxsize:
            # Find all pixels for label
            where = numpy.where(labels == i)
            labels[where] = 0
    # Get set of unique labels
    ulabels = numpy.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[numpy.where(labels == j)] = i
    return labels

labels = filter_by_size(labels, 500, 10000000)

plt.imshow(labels == 25)
plt.show()

kernel = numpy.zeros((9, 9), numpy.float32)
kernel += scipy.stats.norm.pdf(numpy.linspace(-2, 2, 9))[:, na]
kernel += scipy.stats.norm.pdf(numpy.linspace(-2, 2, 9))[na, :]

plt.imshow(kernel)
plt.show()

blurred = scipy.ndimage.convolve(img2, kernel)
plt.imshow(blurred)
plt.show()

img2 = img2 / numpy.amax(img2)
blurred = blurred / numpy.amax(blurred)

blurred2 = numpy.concatenate((img2[:, :, na], blurred[:, :, na]), axis=2)
fig = px.imshow(blurred2, facet_col=2)
fig.show()

kernel = numpy.array([[1, 0, -1], [2, 0, -2], [1, 0, -1]], numpy.float32)
plt.imshow(kernel)
plt.show()

filter1 = scipy.ndimage.convolve(img, kernel)
filter2 = scipy.ndimage.convolve(img, kernel.T)
filter3 = (filter1**2 +filter2**2) ** 0.5
filter1 /= numpy.amax(filter1)
filter2 /= numpy.amax(filter2)
filter3 /= numpy.amax(filter3)
img /= numpy.amax(img)
combined = numpy.concatenate((img[:,:,na],
                              filter1[:, :, na],
                              filter2[:, :, na],
                              filter3[:, :, na]), axis=2)
fig = px.imshow(combined, facet_col=2)
fig.show()
    