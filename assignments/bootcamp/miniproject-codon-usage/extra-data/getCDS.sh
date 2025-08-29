#!/bin/bash

# pip install bio
# conda install jq

bio search -l 100 go:0005737 --species human | \
    jq -r '.[].refseq.rna | if type == "array" then .[0] else . end' | \
    bio fetch | \
    bio fasta --type CDS > cytoplasm.fa

# Address nine 'null's that are returned
bio search -l 109 go:0016020 --species human | \
    jq -r '.[].refseq.rna | if type == "array" then .[0] else . end' | \
    grep -v "null" | \
    bio fetch | \
    bio fasta --type CDS > membrane.fa
