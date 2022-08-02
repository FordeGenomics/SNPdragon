#!/usr/bin/env python

import argparse
from sample import *
import pickle

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str, metavar='sample', help="Name of sample")
    parser.add_argument('-c', type=str, metavar='sample.pileup',  help="Coverage pileup file")
    parser.add_argument('-b', type=str, metavar='sample.bcf', help="Variant bcf file")
    parser.add_argument('--remove_clusters', type=str, choices=('True', 'true', 'False', 'false'), default=False, help='Remove SNPs in clusters (3 SNPs in 10bp window)')
    parser.add_argument('--remove_cliffs', type=str, choices=('True','true','False','false'), default=False, help="Perform cliff searching")
    return parser.parse_args()

def main():
    args = parse_args()
    print("Loading",args.s)

    remove_clusters = args.remove_clusters.lower() == 'true'
    remove_cliffs = args.remove_cliffs.lower() == 'true'


    pileup="{}.pileup".format(args.s)
    bcf="{}.bcf".format(args.s)
    sample = Sample(name=args.s, cov_file=args.c, vcf_file=args.b, filt_clust=remove_clusters, filt_cliffs=remove_cliffs)
    outfile = open("{}.pkl".format(args.s), 'wb')
    pickle.dump(sample, outfile)

if __name__ == '__main__':
    main()
