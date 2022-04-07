#!/usr/bin/env python3

import argparse
from sample import *
import pickle

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str, metavar='sample', help="Name of sample")
    parser.add_argument('-c', type=str, metavar='sample.pileup',  help="Coverage pileup file")
    parser.add_argument('-b', type=str, metavar='sample.bcf', help="Variant bcf file")
    parser.add_argument('-f', type=float, default=0.5, metavar='0.5', help='Fraction of support for snp allele')
    parser.add_argument('-r', action='store_true', default=False, help='Use inbuilt recombination removal')
    parser.add_argument('--keep_cliffs', action='store_false', default=True, help="Turn off cliff searching")
    return parser.parse_args()

def main():
    args = parse_args()
    print("Loading",args.s)

    pileup="{}.pileup".format(args.s)
    bcf="{}.bcf".format(args.s)
    sample = Sample(name=args.s, frac=args.f, cov_file=args.c, vcf_file=args.b, filt_clust=args.r, filt_cliffs=args.keep_cliffs)
    outfile = open("{}.pkl".format(args.s), 'wb')
    pickle.dump(sample, outfile)

if __name__ == '__main__':
    main()
