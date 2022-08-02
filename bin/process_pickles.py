#!/usr/bin/env python

"""
Example: python process_pickles.py -s test/samples.list -o test -t 3 -f 0.5 -r

s : file of sample names
l : list of samples names in format "[sample1, sample2]" # command used for nextflow version
o : output directory
t : number of threads/processes
m : bedfile of positions to mask (e.g. recombination positions output by Gubbins)
r : use inbuilt recombination removal (SNP clusters within sliding window)
 """

import argparse
import os
import subprocess
import re
import processes_nf as pr
import datetime # for testing speed of processes
import numpy as np
import multiprocessing
from sample import *
import pickle
# import tracemalloc

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, required=True, metavar='reference.fa', help="Reference fasta file")
    parser.add_argument('-s', type=str, required=True, metavar='Sample_name', help="Sample name to produce fasta alignment files for")
    parser.add_argument('-i', type=str, metavar='samples_list.txt', help="File of sample names")
    parser.add_argument('-l', type=str, metavar='"[Sample1, Sample2]"', help="list of sample names in the collection utilised by nextflow version")
    # parser.add_argument('-o', type=str, required=True, metavar='path/to/output/dir', help="Output directory")
    parser.add_argument('-t', type=int, default=1,   metavar='1', help='Number of threads/cpus')
    parser.add_argument('-m', type=str, metavar='recomb_pos.gff', help="GFF or Bed file of ranges e.g. recombination and/or phage to mask in the alignment")
    parser.add_argument('-e', type=bool, default=False, help='Exclude reference sequence from output alignment files')
    parser.add_argument('-c', type=int, default=70, metavar='70', help="Minimum percent contig breadth of coverage of sample to include in output")
    parser.add_argument('--remove_clusters', type=str, choices=('True', 'true', 'False', 'false'), default=False, help='Remove SNPs in clusters (3 SNPs in 10bp window)')
    parser.add_argument('--remove_cliffs', type=str, choices=('True','true','False','false'), default=False, help="Perform cliff searching")

    return parser.parse_args()

def read_input(infile):
    samples = []
    with open(infile) as input:
        for line in input:
            samples.append(line.strip())
    return samples

def read_input_file(infile):
    samples = []
    with open(infile) as input:
        for line in input:
            sample = line.strip()
            samples.append(sample)
    return(samples)

def read_input_list(inlist):
    samples = []
    for sample in inlist.strip('[').strip(']').split(', '):
        samples.append(sample)
    return(samples)

def load_pickled_samples(sample):
    print("Unpickling",sample,datetime.datetime.now())
    pickle_file = open("{}.pkl".format(sample), 'rb')
    return pickle.load(pickle_file)

def main():
    print("Commencing processing results:", datetime.datetime.now())
    # tracemalloc.start()
    # current, peak = tracemalloc.get_traced_memory()
    # print(f"Memory tracing start {current / 1024**2}MB; Peak was {peak / 1024**2}MB")
    args = parse_args()
    # make sure dirs have '/' at the end
    # if not args.o.endswith('/'):
    #     args.o = args.o+'/'

    print("Loading results:", datetime.datetime.now())

    if args.i:
        samples = read_input_file(args.i)
        pass
    elif args.l:
        samples = read_input_list(args.l)
    else:
        exit()

    # # Multiprocessing method for creating Sample objects (NOT USED)
    # p = multiprocessing.Pool(processes=min(args.t, len(samples))) # don't request more than necessary
    # loaded_samples = p.map(load_pickled_samples, samples)

    # Commented out 23/07/21
    # non-multiprocessing method
    # loaded_samples = []
    # for sample in samples:
    #     loaded_samples.append(load_pickled_samples(sample))

    # results = pr.Results(loaded_samples)    

    # current, peak = tracemalloc.get_traced_memory()
    # print(f"Memory results loaded {current / 1024**2}MB; Peak was {peak / 1024**2}MB")

    # for sample_obj in loaded_samples:
    #     if sample_obj.name == args.s:
    #         sample = sample_obj

    results = pr.Results(min_cov=args.c)
    for pickled_sample in samples:
        sample_obj = load_pickled_samples(pickled_sample)
        results.addSample(sample_obj)
        if sample_obj.name == args.s:
            sample = sample_obj
        

    print("Building Matrix object:", datetime.datetime.now())
    remove_clusters = args.remove_clusters.lower() == 'true'
    remove_cliffs = args.remove_cliffs.lower() == 'true'

    matrix = pr.Mtx(args.f, sample, results, threads = args.t, recomb = remove_clusters, cliff = remove_cliffs, exclude = args.e, maskFile = args.m)

    # current, peak = tracemalloc.get_traced_memory()
    # print(f"Memory matrix built {current / 1024**2}MB; Peak was {peak / 1024**2}MB")

    # print("Writing full fastas files:", datetime.datetime.now())
    # matrix.writeFastaAlns()

    # print("Writing fuzzy SNP alignments:", datetime.datetime.now())
    # matrix.writeFullSnps()

    # print("Writing core SNP alignments:", datetime.datetime.now())
    # matrix.writeCoreSnps()

    # current, peak = tracemalloc.get_traced_memory()
    # print(f"Memory core SNP alignments saved {current / 1024**2}MB; Peak was {peak / 1024**2}MB")

    # # vcf_files = []
    # # for sample in samples:
    # #     vcf_files.append("{}_core.vcf.gz".format(sample))
    # #     command = ['bcftools','sort','-O','z','-o',"{}_core.vcf.gz".format(sample),"{}_core.vcf".format(sample)]
    # #     sorting = subprocess.Popen(command)
    # #     sorting.communicate()
    # #     command = ['bcftools','index',"{}_core.vcf.gz".format(sample)]
    # #     index = subprocess.Popen(command)
    # #     index.communicate()
    # #     rm_old = subprocess.Popen(['rm','-f',"{}_core.vcf".format(sample)])
    # #     rm_old.communicate()
    # #
    # # # current, peak = tracemalloc.get_traced_memory()
    # # # print(f"Memory vcf files sorted {current / 1024**2}MB; Peak was {peak / 1024**2}MB")
    # #
    # # # tracemalloc.reset_peak()
    # # # Merge core vcf files to single multi-sample vcf
    # # command = ['bcftools', 'merge', '-0', '-o', 'multisample_core.vcf']
    # # command.extend(vcf_files)
    # # merge = subprocess.Popen(command)
    # # merge.communicate()
    # #
    # # current, peak = tracemalloc.get_traced_memory()
    # # print(f"Memory vcf files merged {current / 1024**2}MB; Peak was {peak / 1024**2}MB")

if __name__ == '__main__':
    main()
