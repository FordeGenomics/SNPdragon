#!/usr/bin/env python3

import argparse
import os
import subprocess
import datetime
import pandas as pd
# import tracemalloc

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, required=True, metavar='reference.fa', help="Reference fasta file")
    parser.add_argument('-i', type=str, metavar='samples_list.txt', help="File of sample names")
    parser.add_argument('-l', type=str, metavar='"[Sample1, Sample2]"', help="list of sample names in the collection utilised by nextflow version")
    parser.add_argument('-o', type=str, required=True, metavar='path/to/output/dir', help="Output directory")
    parser.add_argument('-e', action='store_true', default=False, help='Exclude reference sequence from output alignment files')

    return parser.parse_args()

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

def calc_snp_mtx(samples):
    first = True
    cols = ['Contig']
    coles = cols.extend(samples)
    for sample in samples:
        filename = '{}_snpmatrix.csv'.format(sample)
        if os.path.exists(filename):
            if first:
                first = False
                data = pd.read_csv(filename, index_col=1)
                mtx = pd.DataFrame(index=list(data.index), columns=cols)
                mtx.update(data)
            else:
                data = pd.read_csv(filename, index_col=1)
                mtx.update(data)
    try:
        mtx.insert(1, "POS", data.index)
        output = open("snp_matrix.csv", "w")
        mtx.to_csv(output, index=False)
    except UnboundLocalError: # mtx never initiated, no samples included
        pass

def main():
    print("Concatenating results", datetime.datetime.now())

    args = parse_args()
    # make sure dirs have '/' at the end
    if not args.o.endswith('/'):
        args.o = args.o+'/'

    if args.i:
        samples = read_input_file(args.i)
    elif args.l:
        samples = read_input_list(args.l)
    else:
        exit()

    cat_core_commands = ['cat']
    cat_full_commands = ['cat']
    cat_stats_commands = ['cat']
    cat_full_aln_commands = ['cat']
    cat_exclude_commands = ['cat']

    calc_snp_mtx(samples)
    
    included = 0
    excluded = 0
    print(samples)
    for sample in sorted(samples):
        if os.path.exists('{}_core_snp.fasta'.format(sample)):
            included += 1
            cat_core_commands.append('{}_core_snp.fasta'.format(sample))
            cat_full_commands.append('{}_full_snp.fasta'.format(sample))
            cat_stats_commands.append('{}_core_stats.tsv'.format(sample))
            cat_full_aln_commands.append('{}_full_aln.fasta'.format(sample))
        if os.path.exists('{}_excluded.tsv'.format(sample)):
            excluded += 1
            cat_exclude_commands.append('{}_excluded.tsv'.format(sample))

    
    print("Included:",included)
    print("Excluded:",excluded)

    if included > 0:
        core_file = open("core_snp.fasta", 'w')
        cat_core = subprocess.Popen(cat_core_commands, stdout=core_file)
        cat_core.communicate()
        core_file.close()

        full_file = open("full_snp.fasta", 'w')
        cat_full = subprocess.Popen(cat_full_commands, stdout=full_file)
        cat_full.communicate()
        full_file.close()

        # core_aln_file = open("core_aln.fasta", 'w')
        # cat_core_aln = subprocess.Popen(cat_core_aln_commands, stdout=core_aln_file)
        # cat_core_aln.communicate()
        # core_aln_file.close()

        full_aln_file = open("full_aln.fasta", 'w')
        cat_full_aln = subprocess.Popen(cat_full_aln_commands, stdout=full_aln_file)
        cat_full_aln.communicate()
        full_aln_file.close()

        stats_file = open("core_stats.tsv", 'w')
        stats_file.write("#ID\tCONTIG\tCORE_SIZE\tCOVERAGE%\n")
        stats_file.flush()
        cat_stats = subprocess.Popen(cat_stats_commands, stdout=stats_file)
        cat_stats.communicate()
        stats_file.close()

        # get snp distance matrix
        print("Calc core dist", datetime.datetime.now())
        dist_file = open("snp_dist.csv", 'w')
        dist_commands = ['snp-dists', '-b', '-c', '-q', 'core_snp.fasta']
        dist = subprocess.Popen(dist_commands, stdout=dist_file)
        dist.communicate()
        dist_file.close()
        print("Finished", datetime.datetime.now())
    
    if excluded > 0:
        exclude_file = open("excluded.tsv", 'w')
        exclude_file.write("#ID\tCONTIG\tCOVERAGE%\n")
        exclude_file.flush()
        cat_exclude = subprocess.Popen(cat_exclude_commands, stdout=exclude_file)
        cat_exclude.communicate()
        exclude_file.close()

if __name__ == '__main__':
    main()
