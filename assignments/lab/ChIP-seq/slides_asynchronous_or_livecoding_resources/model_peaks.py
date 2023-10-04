#!/usr/bin/env python

import sys

import numpy
import matplotlib.pyplot as plt

def main():
    # Load bedgraph datasets
    forward_fname, reverse_fname, out_fname = sys.argv[1:4]

    # Define our variables and target region
    target = "chr2R"
    chromlen = 25286936
    binsize = 200
    number_of_peaks = 200

    # Load in the bedgraph data
    forward = load_bedgraph(forward_fname, target, 0, chromlen)
    reverse = load_bedgraph(reverse_fname, target, 0, chromlen)

    # Combine forward and reverse tags
    combined = forward + reverse

    # Bin tag density using a sliding window
    scores = bin_array(combined, binsize)

    # Identify the top N peak positions
    peaks = find_peaks(scores, number_of_peaks)

    # Find combined tag densities over the selected peaks by strand
    reverse_curve = find_profile(reverse, peaks, binsize * 4)
    forward_curve = find_profile(forward, peaks, binsize * 4)

    # Find best shift to match forward and reverse peaks
    correlations = find_correlations(reverse_curve, forward_curve)

    # Plot results
    fig, ax = plt.subplots(2, 1, figsize=(5, 10))

    # Plot unshifted strand profiles
    X = numpy.arange(-reverse_curve.shape[0] // 2, reverse_curve.shape[0] // 2)
    ax[0].plot(X, reverse_curve, color='blue', label='reverse')
    ax[0].plot(X, forward_curve, color='red', label='foward')
    ax[0].set_xlabel("Distance from combined peak (bp)")
    ax[0].set_ylabel("Mean tag count")
    ax[0].set_title("Peak offset by strand")

    # Plot correlation curve
    ax[1].plot(numpy.arange(correlations.shape[0]), correlations, color='black')
    ax[1].set_xlabel("Fragment size")
    ax[1].set_ylabel("Correlation")
    ax[1].set_title('Correlation by offset')

    # Mark correlation peak
    best = numpy.argmax(correlations)
    bestval = numpy.amax(correlations)
    ax[1].plot([best, best], [0, bestval], color='black')
    plt.tight_layout()
    plt.savefig(out_fname)
    plt.close()

    # Report ideal shift for downstream analysis
    print("Best offset was {}".format(best))


def load_bedgraph(fname, target, chromstart, chromend):
    # Create array to hold tag counts
    coverage = numpy.zeros(chromend - chromstart, int)

    #Read the file in line by line
    for line in open(fname):
        # Break the line into individual fields
        chrom, start, end, score = line.rstrip().split('\t')
        # Check if the data fall in our target region
        if chrom != target:
            continue
        start = int(start)
        end = int(end)
        if start < chromstart or end >= chromend:
            continue
        # Add tags to our array
        coverage[start-chromstart:end-chromend] = int(score)
    return coverage

def bin_array(data, binsize):
    # Create array to hold scores
    binned = numpy.zeros(data.shape[0], data.dtype)

    # For each position in the window, add to the score array
    for i in range(binsize):
        binned[i:data.shape[0] - binsize + i] += data[binsize//2:-binsize//2]
    return binned

def find_peaks(data, target):
    # Find the ordering of peak scores are reverse (high to low)
    order = numpy.argsort(data)[::-1]
    # Return only the desired number of peaks
    return order[:target]

def find_profile(data, peaks, binsize):
    # Create array to hold the combined profile
    results = numpy.zeros(binsize, float)
    for i in range(peaks.shape[0]):
        # For each peak, get the tag density +/- half the binsize
        results += data[peaks[i]-binsize//2:peaks[i]+binsize//2]
    # Return average tag density
    return results / peaks.shape[0]

def find_correlations(reverse, forward):
    # Since we made the profiles 4 times bigger than our estimated fragment size
    # we will search a size from 0 to 2 times our estimate
    width = reverse.shape[0] // 2
    # Create an array to hold correlations
    corrs = numpy.zeros(width, float)
    for i in range(width):
        # For each shift, find correlation of overlapping regions 
        corrs[i] = numpy.corrcoef(reverse[i:], forward[:forward.shape[0]-i])[0, 1]
    return corrs

# Make sure that we can import things from this script without it runnning
if __name__ == "__main__":
    main()
