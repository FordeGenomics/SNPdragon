#!/usr/bin/env python

import sys

fasta_index_file = sys.argv[1]
if not fasta_index_file.endswith(".fai"):
    fasta_index_file = fasta_index_file + ".fai"

fasta_index_file = open(fasta_index_file)
region_size = int(sys.argv[2])

for line in fasta_index_file:
    fields = line.strip().split("\t")
    chrom_name = fields[0]
    chrom_length = int(fields[1])
    region_start = 0
    while region_start < chrom_length:
        start = region_start
        end = region_start + region_size -1
        if end > chrom_length:
            end = chrom_length
        print(chrom_name + ":" + str(region_start) + "-" + str(end))
        region_start += region_size
