#!/bin/bash

#USAGE: bash exercise4.sh

#Goal: report any genes that uniquely intersect with H3K27ac but never intersect with H3K9me3 within naive B cells.

bedtools intersect -a H3K27ac.naive_b_cell.GRCh38.bedgraph -b genes.bed -wb | cut -f 8 | sort | uniq > genes_intersecting_H3K27ac_b_cell.txt

bedtools intersect -a H3K9me3.naive_b_cell.GRCh38.bedgraph -b genes.bed -wb | cut -f 8 | sort | uniq >
genes_intersecting_H3K9me3_b_cell.txt

grep -v -f genes_intersecting_H3K27ac_b_cell.txt genes_intersecting_H3K9me3_b_cell.txt
