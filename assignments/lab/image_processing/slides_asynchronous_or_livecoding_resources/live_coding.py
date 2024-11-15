#!/usr/bin/env python

import numpy
import scipy
import matplotlib.pyplot as plt
import imageio

na = numpy.newaxis

# Load image into a numpy 2D array and rescale
img = imageio.v3.imread("illum_DAPI.tif").astype(numpy.float32)
img -= numpy.amin(img)
img /= numpy.amax(img) * 1.05

img1 = imageio.v3.imread("illum_RNAcytoNuc.tif").astype(numpy.float32)
img1 -= numpy.amin(img1)
img1 /= numpy.amax(img1) * 1.05

img2 = imageio.v3.imread("illum_Mito.tif").astype(numpy.float32)
img2 -= numpy.amin(img2)
img2 /= numpy.amax(img2) * 1.05

# Check image size and data type
print(img.shape, img.dtype)

# Display image
plt.imshow(img, vmin=0, vmax=2**16)
plt.show()

# Look at data range
plt.hist(img)
print(numpy.amin(img), numpy.amax(img))

# Display image with automatic value scaling
plt.imshow(img)
plt.show()

# Load 3 channels into single image array
rgbimg = numpy.zeros((img.shape[0], img.shape[1], 3), numpy.uint16)
for i, name in enumerate(['DAPI', 'RNAcytoNuc', 'Mito']):
    rgbimg[:, :, i] = imageio.v3.imread(f"illum_{name}.tif")

# Display image
plt.imshow(rbgimg)

# Oops, got a warning.
# Convert values from uint16 to uint8
rgbimg = (rgbimg // 2**8).astype(numpy.uint8)

# Display image
plt.imshow(rbgimg)
plt.show()

# Too dim because there is no automatic color scaling
# Maximize range for each color
for i in range(3):
    minv = numpy.amin(rgbimg[:, :, i])
    maxv = numpy.amax(rgbimg[:, :, i])
    newimg = (rgbimg[:, :, i] - minv) / (maxv - minv)
    newimg = numpy.round(newimg * 240).astype(numpy.uint8)
    rgbimg[:, :, i] = newimg

# Display image
fig = px.imshow(rgbimg, facet_col=2)
fig.show()

# Let's look at creating a masking image
mask = img > 0.1
plt.imshow(mask)
plt.show()

# Let's use a function to find nuclei bodies
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
# Since the first label is 1 and the background is 0, let's adjust the background for more contrast
label_copy = numpy.copy(labels)
label_copy[numpy.where(label_copy == 0)] -= 50
plt.imshow(label_copy)
plt.show()

# # Filter cells
sizes = numpy.bincount(labels)
print(sizes)
for i in range(sizes.shape[0]):
    if sizes[i] < 100:
        labels[numpy.where(labels == i)] = 0

def filter_by_size(labels, minsize, maxsize):
    # Find label sizes
    sizes = numpy.bincount(labels)
    # Iterate through labels, skipping background
    for i in range(1, sizes.shape[0]):
        # If the number of pixels falls outsize the cutoff range, relabel as background
        if sizes[i] < minsize or sizes > maxsize:
            # Find all pixels for label
            where = numpy.where(labels == i)
            labels[where] = 0
    # Get set of unique labels
    ulabels = numpy.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[numpy.where(labels == j)] = i
    return labels

# Let's look at the labels after filtering
label_copy = numpy.copy(labels)
label_copy[numpy.where(label_copy == 0)] -= 50
plt.imshow(label_copy)
plt.show()

# What if we want to select a single marked nucleus?
marked = numpy.copy(mask).astype(numpy.int32)
where = numpy.where(label == 50)
marked[where] = 2

# And if we want information about that nucleus in another channel?
minv = numpy.amin(img1[where])
maxv = numpy.amax(img1[where])
print(f"RNA nuclear signal ranges from {minv} to {maxv}")

# Now let's see what a kernel is and what it does
kernel = numpy.zeros((9, 9), numpy.float32)
# Add two normal curves, one across rows, one across columns
kernel += scipy.stats.norm.pdf(numpy.linspace(-2, 2, 9))[:, na]
kernel *= scipy.stats.norm.pdf(numpy.linspace(-2, 2, 9))[na, :]

plt.imshow(kernel)
plt.show()

# Let's see what happens when we apply the kernel
blurred = scipy.ndimage.convolve(img2, kernel)

# Let's try a new way of seeing the data
import plotly.express as px
import plotly

# We need to combine the images
combined = numpy.concatenate((img2[:, :, na], blurred[:, :, na]), axis=2)

# Now we can view them with plotly
fig = px.imshow(combined, facet_col=2)
fig.show()

# Let's look at another kernel
kernel = numpy.array([[1, 0, -1], [2, 0, -2], [1, 0, -1]])
plt.imshow(kernel, vmin=-2, vmax=2)
plt.show()

# Apply kernel in vertical and horizontal directions
filtered_x = scipy.ndimage.convolve(img, kernel)
filtered_y = scipy.ndimage.convolve(img, kernel.T)
filtered = (filtered_x**2 + filtered_y**2) ** 0.5

# Now let's combine the images so we do a fancy display
fullfig = numpy.concatenate((img[:,:,na], filtered_x[:,:,na], filtered_y[:,:,na], filtered[:,:,na]), axis=2)

# And rescale each layer independently
fullfig -= numpy.amin(fullfig.reshape(-1, 4), axis=0)
fullfig /= numpy.amax(fullfig.reshape(-1, 4), axis=0)

# And finally plot them in a linked and interactive plot
fig = px.imshow(fullfig, facet_col=2)
fig.show()

