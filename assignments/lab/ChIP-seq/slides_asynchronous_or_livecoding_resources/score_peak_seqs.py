#!/usr/bin/env python

import sys

import numpy
import matplotlib.pyplot as plt


def main():
    # Load input files
    
    # Define a pseudocount to avoid taking the log of zero
    pseudocount = 1.0
    
    # Load the position frequency matrix (pfm)

    # Load the peak coordinates from a bed file. Use the size of the pfm for width

    # Load the sequence for our chromosome

    # For each peak, find the best score and best scoring sequence

    # Plot histogram of scores, original pfm, and the mean pfm from the peak sequences
    # Note that because it's a 2D matrix, you can just use imshow to plot the pfms.

def load_bed(fname, width):
    bed = []
    for line in open(fname):
        line = line.rstrip().split('\t')
        start = int(line[1])
        end = int(line[2])
        if end - start < width:
            continue
        bed.append((start, end))
    bed = numpy.array(bed, int)
    return bed

def load_fasta(fname):
    seq = []
    for line in open(fname):
        if line.startswith(">"):
            continue
        seq.append(line.rstrip())
    seq = "".join(seq).upper()
    return seq

def load_pfm(fname, pseudocount):
    pfm = []
    for line in open(fname):
        if line.startswith(">"):
            continue
        pfm.append(line.rstrip().split('\t'))
    pfm = numpy.array(pfm, float) + pseudocount
    pfm /= numpy.sum(pfm, axis=0, keepdims=True)
    return pfm

def seq_to_onehot(seq):
    onehot = numpy.zeros((4, len(seq)), int)
    index = {"A": 0, "C": 1, "G": 2, "T": 3}
    for i in range(len(seq)):
        onehot[index[seq[i]], i] = 1
    return onehot

def score_seqs(peaks, seq, pfm, pseudocount):
    # Convert the pfm into a pwm by log10 transformation
    
    # Get the motif size (K) from the pwm or pfm shape[1]

    # Calculate what score of a sequence would be if all bases were equally likely
    # Think a pfm full of 0.25 values

    # Create a way to store your results

    #Iterate through each peak
        # Extract the sequence corresponding to the particular peak

        # skip peak sequences with Ns

        # For each peak sequence, convert it into onehot encoding

        # Create a temporary array to hold the scores. Think how long the array
        # needs to be since each score requires a K-mer worth of bases
        # Also remember that the motif could be in either orientation
        
        # Iterate through each possible k-mer in the peak sequence

            # For each K-mer, score the onehot sequence by multiplying the pwm
            # and onehot sequence and adding up the total

            # Score the reverse-compliment the same way. Note that to reverse-
            # complement a onehot encoded sequence, you just need to reverse the
            # order of both axes

        # Check which strand had the best score 
        # Record the best score and add the best onehot sequence to your mean pfm

    # Subtract the background from your scores (since these are log-transformed)
    # this represents how many times more likely the motif is over background at
    # the sequence

    # Add the pseudocounts

    #Convert your mean pfm from counts to percentages (divide by the sum along columns)

    # Return the scores and pfm obtained from the peak sequences



if __name__ == "__main__":
    main()