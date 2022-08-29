#!/bin/bash

#USAGE: bash exercise3.sh input_VCF

awk '/^#/{next} {print $1,$2-1, $2}' $1 > variants.bed
sort -k1,1 -k2,2n genes.bed > genes.sorted.bed
bedtools closest -a variants.bed -b genes.sorted.bed | head
