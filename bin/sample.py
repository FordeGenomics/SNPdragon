#!/usr/bin/env python3

import csv
import numpy as np
import pandas as pd
from collections import OrderedDict
# import pyranges as pyr
from pysam import VariantFile
import scipy.sparse as sp
from statistics import mean
import multiprocessing
import os
import sys
import datetime
from dna_dicts import *
from json import JSONEncoder

# import tracemalloc

import json

# Globals
MIN_DEPTH = 10
SNP_CLUST = 3
WINDOW = 10
ALT_RATIO = 0.75
CLIFF_WINDOW = 10
CLIFF_SLOPE = 3
SEARCH_WINDOW = 10
R_SQUARED_THRESHOLD = 0.7
MIN_MQM = 40
READ_BALANCE = 0.05

class Position:
    def __init__(self, contig, pos):
        self.contig = contig
        self.pos = int(pos)

    def __str__(self):
        return "{}:{}".format(self.contig, self.pos)

    def __repr__(self):
        return "{}:{}".format(self.contig, self.pos)

    def __eq__(self, other):
        if isinstance(other, Position):
            return self.contig == other.contig and self.pos == other.pos
        elif isinstance(other, int):
            return self.pos == other
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, Position):
            return self.contig != other.contig or self.pos != other.pos
        elif isinstance(other, int):
            return self.pos != other
        else:
            return NotImplemented

    def __lt__(self, other):
        if isinstance(other, Position):
            if self.contig == other.contig:
                return self.pos < other.pos
            else:
                return self.contig < other.contig
        elif isinstance(other, int):
            return self.pos < other
        else:
            return NotImplemented

    def __le__(self, other):
        if self.__eq__(self, other):
            return True
        else:
            return self.__lt__(self, other)

    def __gt__(self, other):
        if isinstance(other, Position):
            if self.contig == other.contig:
                return self.pos > other.pos
            else:
                return self.contig > other.contig
        elif isinstance(other, int):
            return self.pos > other
        else:
            return NotImplemented

    def __ge__(self, other):
        if self.__eq__(self, other):
            return True
        else:
            return self.__gt__(self, other)

    def __hash__(self):
        return hash("{}:{}".format(self.contig,self.pos))

    def __int__(self):
        """ Allow casting to type int """
        return self.pos

    def __float__(self):
        return float(self.pos)

    def __add__(self, other):
        return self.pos + other.pos

    def __sub__(self, other):
        return self.pos - other.pos

    def sub(self, n):
        return self.pos - n

    def add(self, n):
        return self.pos + n

class Coverage(Position):
    def __init__(self, contig, pos, depth):
        self.depth = depth
        super().__init__(contig, pos)

class SNP(Position):
    """ Variant class storing details of a single variant called.
    Greater than, equal to, less than with int class also implemented to allow
    for finding matching positions of coverage
    @TODO: Option of implementing a Variant SNP and COV classes instead of allowing
    multiple class comparisons
    params:
        contig (str): contigosome/reference name variant was detected in
        pos (int): start position
        alt (str): alternate allele (if SNP or small indel, otherwise '*') """

    def __init__(self, contig, pos, alt):
        self.alt = alt
        super().__init__(contig, pos)
    
    def float(self):
        return float("{}.{}".format(self.pos, self.alt))

    def __str__(self):
        return "{}:{}:{}".format(self.contig, self.pos, self.alt)

    def __repr__(self):
        return "{}:{}:{}".format(self.contig, self.pos, self.alt)

    def __eq__(self, other):
        if isinstance(other, SNP):
            return self.contig == other.contig and self.pos == other.pos and self.alt == other.alt
        elif isinstance(other, str):
            return self.__str__() == other
        else:
            return Position.__eq__(self, other)

    def __ne__(self, other):
        if isinstance(other, SNP):
            return self.contig != other.contig or self.pos != other.pos
        elif isinstance(other, str):
            return self.__str__() != other
        else:
            return Position.__ne__(self, other)

    def __lt__(self, other):
        if isinstance(other, SNP):
            if self.contig == other.contig:
                if self.pos == other.pos:
                    return self.alt < other.alt
                else:
                    return self.pos < other.pos
            else:
                return self.contig < other.contig
        elif isinstance(other, str):
            contig, pos, alt = other.split(":")
            if self.contig == contig:
                if self.pos == int(pos):
                    return self.alt < alt
                else:
                    return self.pos < int(pos)
            else:
                return self.contig < contig
        else:
            return Position.__lt__(self, other)

    def __le__(self, other):
        if self.__eq__(self, other):
            return True
        else:
            return self.__lt__(self, other)

    def __gt__(self, other):
        if isinstance(other, SNP):
            if self.contig == other.contig:
                if self.pos == other.pos:
                    return self.alt > other.alt
                else:
                    return self.pos > other.pos
            else:
                return self.contig > other.contig
        elif isinstance(other, str):
            contig, pos, alt = other.split(":")
            if self.contig == contig:
                if self.pos == int(pos):
                    return self.alt > alt
                else:
                    return self.pos > int(pos)
            else:
                return self.contig > contig
        else:
            return Position.__gt__(self, other)

    def __ge__(self, other):
        if self.__eq__(self, other):
            return True
        else:
            return self.__gt__(self, other)

    def __hash__(self):
        return hash("{}:{}:{}".format(self.contig,self.pos,self.alt))

    def getShortRepr(self):
        return "{}:{}".format(self.contig, self.pos)

class Sample:
    """ Name, coverage csv from bedfile and snp results for each sample
    Attr:
        name (str): Sample name
        snps (dict): Reference allele at each SNP position {position:reference}
        snpSeries (Pandas Series): Filtered/high quality final SNPs calls, colname = position, value = alt allele
        ufSeries (Pandas Series): Unfiltered SNP calls (will be same as snpSeries if only 1 vcf file passed)
        covSeries (Pandas Series): Low coverage positions below depth threshold,
                                    colnames = position, value = N
        allCov (Pandas Dataframe): [contig, Pos, Depth]
        cliffs (list): List of positions occuring in cliffs
    """

    def __init__(self, name, frac, cov_file, vcf_file, outdir="", filt_clust=False, filt_cliffs=True):
        """ Name, coverage csv from bedfile and snp results for each sample
        Params:
            name (str): Sample name
            cov_file (str): File for per base coverage output by bedtools -d
            vcf_filtered (str): VCF file of all SNP calls
        """
        self.name = name
        self.outdir = outdir
        self.frac = frac # Default 0.5 passed in pickle_samples.py TO BE REMOVED
        self.filter_clust = filt_clust
        self.filter_cliffs = filt_cliffs
        self.contigs = set()
        self.header = None
        self.variants = OrderedDict() # stores all variantFile records, cleaned
        self.refsUnfilt = {} # stores reference allele of all variant positions key = str(CONTIG:POS)
        self.snps = {} # cleaned
        self.snpsUnfilt = {} # cleaned
        self.snpSeries = {}
        self.ufSeries = {}
        self.covSeries = {}
        self.allCov = {} # cleaned
        self.cliffs = {}
        self.mix = {} # set of low freqency variants below fraction of support
        self.clusters = {}
        self.coverage = {}

        self.readPileup(cov_file)
        self.readVCF(vcf_file)

        if self.filter_clust:
            self.findClusters()

        if self.filter_cliffs:
            self.findCliffs()

        self.buildSeries()
        self.writeVcf()

        # cleanup
        del self.allCov
        del self.snps
        del self.snpsUnfilt
        del self.variants


    def readVCF(self, file):
        """ Read in VCF file of SNPs
        TODO: Read in variant positions with multiple alleles to save best quality
        TODO: Handle multiple contigosome/contig reference """
        print("Reading in vcf for", self.name, datetime.datetime.now())
        vcf_in = VariantFile(file, mode='r')
        vcf_in.header.add_meta('FORMAT', items=[\
            ('ID',"FT"), 
            ('Number',"."), 
            ('Type',"String"), 
            ('Description',"Filtering for inclusion in matrix. PASS:AO/DP>={} AND DP>={} AND MQM>{}.".format(ALT_RATIO, MIN_DEPTH, MIN_MQM))
            ])

        self.header = str(vcf_in.header)

        for rec in vcf_in.fetch():
            contig = rec.contig
            self.contigs.add(contig)
            if not self.variants.get(contig):
                self.variants[contig] = OrderedDict()
                self.refsUnfilt[contig] = {}
                self.snps[contig] = {}
                self.snpsUnfilt[contig] = {}
                self.mix[contig] = set()

            multi = len(rec.alts) > 1
            prev_allele = ""
            for allele in range(0,len(rec.alts)): # if mnp or complex (should be none if processed by vcfallelicprimitives)
                if len(rec.alts) > 1: # if multiple alternate calls depth becomes a tuple
                    depth = rec.info['DP'][allele]
                else:
                    depth = rec.info['DP']
                alt_ratio = rec.info['AO'][allele]/depth >= ALT_RATIO
                min_depth = depth >= MIN_DEPTH
                min_mqm = rec.info['MQM'][allele] >= MIN_MQM
                read_balance = \
                    rec.info['SAF'][allele]/(rec.info['SAF'][allele] + rec.info['SAR'][allele]) >= READ_BALANCE \
                        and rec.info['SAR'][allele]/(rec.info['SAF'][allele] + rec.info['SAR'][allele]) >= READ_BALANCE
                qual = alt_ratio and min_depth and min_mqm and read_balance
                if len(rec.ref) == len(rec.alts[allele]): # snp, mnp or complex of equal length to reference
                    for e, r in enumerate(rec.ref):
                        filter = []
                        curr_allele = DNA_CODE[rec.alts[allele][e].upper()]
                        lower_allele = DNA_CODE[rec.alts[allele][e].lower()]
                        var = SNP(contig,rec.pos+e,curr_allele)
                        ref = DNA_CODE[r]
                        # if multi:
                        #     curr_allele = IUPAC[prev_allele+rec.alts[allele][e]]
                        # else:
                        #     curr_allele = rec.alts[allele][e]
                        if ref != curr_allele:
                            self.refsUnfilt[contig][rec.pos] = DNA_CODE[r]
                            self.snpsUnfilt[contig][var] = lower_allele

                            if qual:
                                filter.append("PASS")
                                self.snps[contig][var] = curr_allele
                            else:
                                # set lower than 0.5 (majority rules) for populating matrix to allow for some variation (balance of false positive)
                                if not alt_ratio and rec.info['AO'][allele]/depth >= 0.5: 
                                    filter.append("FAIL_AF")
                                if rec.info['AO'][allele]/depth < 0.5:
                                    filter.append("FAIL_AF0.5")
                                    self.mix[contig].add(rec.pos)
                                    self.snpsUnfilt[contig][var] = DNA_CODE['N']
                                if not min_depth:
                                    filter.append("FAIL_DEPTH")
                                if not min_mqm:
                                    filter.append("FAIL_MQM")
                                if not read_balance:
                                    filter.append("FAIL_RB")
                            if multi:
                                prev_allele = prev_allele+rec.alts[allele][e]

                            rec.samples[0]["FT"] =  ','.join(filter)
                            self.variants[contig][rec.pos+e] = str(rec)

                        
                else: #indel
                    # mask up 2bp either side of indel edges (and across entire del region i.e. len(rec.ref))
                    for n in range(max(rec.pos-1,0),rec.pos+len(rec.ref)+2):
                        self.covSeries[contig].append(n)

        for contig in self.covSeries.keys(): # remove any duplicate pos
            self.covSeries[contig] = sorted(list(set(self.covSeries[contig])))

    def readPileup(self, file):
        """ Read in pileup file of coverage for file generated by samtools depth """
        print("Reading in pileup for", self.name, datetime.datetime.now())
        allCov = pd.read_csv(file, sep="\t", names = ["Contig", "Pos", "Depth"])
        contigs = allCov.Contig.unique().tolist()
        print("Manipulating pileup")
        for contig in contigs:
            self.contigs.add(contig)
            self.allCov[contig] = allCov.loc[allCov["Contig"] == contig][["Depth"]]
            self.allCov[contig].index = allCov.loc[allCov["Contig"] == contig]["Pos"].tolist()
            self.allCov[contig]["Depth"] = self.allCov[contig].loc[:,("Depth")].astype('Int16')

            print("Removing positions above depth threshold", self.name, contig)

            # Remove positions above depth threshold
            self.covSeries[contig] = list(self.allCov[contig].loc[self.allCov[contig]["Depth"] < MIN_DEPTH].index)

            # calculate percent of coverage of genome
            genome_size = len(self.allCov[contig].index)
            low_cov_size = len(self.covSeries[contig])
            self.coverage[contig] = round(100 - low_cov_size/genome_size * 100, 2)

            cov_out = open("{}_{}_low_cov.txt".format(self.name, contig), 'w')
            cov_out.write(",".join([str(n) for n in self.covSeries[contig]]))

    def findClusters(self):
        """ Remove SNPs occuring in clusters. 
        
        Only searches SNPs passed previous filtering ."""
        print("Finding clusters for", self.name, datetime.datetime.now())
        for contig in self.snps.keys():
            clustered = self.clustered(self.snps[contig].keys())
            self.clusters[contig] = set()
            # Update filter in clustered SNPs
            for var in clustered:
                self.variants[contig][var.pos] = self.variants[contig][var.pos].replace("PASS", "FAIL_CLUST")
                self.snps[contig].pop(var)
                self.clusters[contig].add(var.pos)

    def clustered(self, snps):
        """ Returns SNP positions removing clusters of SNP_CLUST occuring within WINDOW in single sample
        Params:
            snps (list): List of snp positions to screen """
        sorted_pos = sorted(snps)
        recomb = set()
        end = len(sorted_pos)
        start = 0
        while start <= end-SNP_CLUST:
            # E.g. SNPs of 3 positions apart in list fall within 10bp difference of each other
            if sorted_pos[start+SNP_CLUST-1] - sorted_pos[start] < WINDOW:
                for n in range(0, SNP_CLUST):
                    recomb.add(sorted_pos[start+n])
            start += 1
        return recomb

    def findCliffs(self):
        """ Search for cliffs in SEARCH_WINDOW either side of a list of positions
        Reports cliffs that are within SEARCH_WINDOW either side of a position if the slope
        is > or < CLIFF_SLOPE and the fit of the line is > R_SQUARED_THRESHOLD

        Only searches SNPs that passed previous filtering

        TODO: Check SNPs occuring in sudden drop-offs of depth and if cliffs are reported at these positions """
        print("Finding cliffs for", self.name, datetime.datetime.now())
        x = np.array([i for i in range(1, CLIFF_WINDOW+1)], dtype=np.float64) # 1-based x coordinates for graph
        total_cliffs = 0
        for contig in self.contigs:
            self.cliffs[contig] = set()
            allCov = self.allCov[contig]
            cov = np.array(list(allCov["Depth"]), dtype=np.float64)
            i = 0
            if self.snps.get(contig) is not None:
                vars = list(self.snps.get(contig).keys())
                for var in vars:
                    start = max(var.pos - SEARCH_WINDOW, 0) # don't go below 0
                    end = min(var.pos + SEARCH_WINDOW, len(cov)) # don't go beyond end of array
                    while start < end-CLIFF_WINDOW: # search sliding window
                        y = cov[start:start+CLIFF_WINDOW]
                        m, b = best_fit_slope_and_intercept(x,y)
                        if m >= CLIFF_SLOPE or m <= -CLIFF_SLOPE: # if slope, get line fit
                            regression_line = [(m*i)+b for i in x]
                            r_squared = coefficient_of_determination(y,regression_line)
                            if r_squared > R_SQUARED_THRESHOLD:
                                self.cliffs[contig].add(var.pos)
                                if self.variants[contig].get(var.pos):
                                    if "PASS" in self.variants[contig][var.pos]:
                                        self.variants[contig][var.pos] = \
                                            self.variants[contig][var.pos].replace("PASS", "FAIL_CLIFF")
                                    else:
                                        self.variants[contig][var.pos] = \
                                            self.variants[contig][var.pos].replace(":FAIL", ":FAIL_CLIFF,FAIL")
                                if self.snps[contig].get(var):
                                    self.snps[contig].pop(var)
                                total_cliffs += 1
                                break
                        start += 1
                    i += 1
            print("Number of cliffs for", self.name ,total_cliffs, datetime.datetime.now())

    def buildSeries(self):
        for contig in self.snps.keys():
            self.snpSeries[contig] = pd.DataFrame([self.snps[contig].values()], index = [\
                self.name], columns=[snp.float() for snp in self.snps[contig].keys()], dtype='uint8')
            self.snpSeries[contig].to_csv("{}_snpSeries.csv".format(self.name)) # REMOVE IN PROD
            ufSeries = pd.DataFrame([self.snpsUnfilt[contig].values()], index = \
                [self.name], columns=[snp.float() for snp in self.snpsUnfilt[contig].keys()], dtype='uint8')
            self.ufSeries[contig] = ufSeries

        for contig in self.covSeries.keys(): # coverage may have contigs with no snps called
            self.covSeries[contig] = np.array(self.covSeries[contig], dtype='uint32')

    def getFiltered(self):
        """ Get low coverage regions and high quality SNPs.
        If a pos is in both low coverage and unfiltered vcf call, low coverage is returned.
        concat([coverage,unfiltered]) in this order to keep coverage sites in removeDupCols """
        series = {}
        for contig in self.contigs:
            if self.snpSeries.get(contig) is None:
                series[contig] = pd.DataFrame([],index=[self.name])
            else:
                series[contig] = self.snpSeries.get(contig)
        return series

    def getUnfiltered(self, contig, ref_series, ref_pos_series):
        """ Return sub-series of unfiltered VCF calls and low coverage from list of SNP positions to populate matrix.
        If a pos is in both low coverage and unfiltered vcf call, it is populated as low coverage """
        # Create low coverage series
        all_pos = ref_series.columns.to_list()
        low_cov_cols = self.getLowCov(all_pos, contig)
        low_cov_vars = pd.DataFrame(np.zeros((1,len(low_cov_cols))), index=[self.name], columns=low_cov_cols, dtype='uint8')

        if self.snpSeries.get(contig) is not None:
            # Prepare required reference series
            ref_pos_series.index = [self.name]

            print("Getting variants for",self.name,datetime.datetime.now())
            
            # Get all passed SNPs
            snpSeries = self.snpSeries[contig][set(all_pos).intersection(set(self.snpSeries[contig].columns))]

            ufSeries = self.ufSeries[contig][set(all_pos).intersection(set(self.ufSeries[contig].columns))]

            # get reference for missing positions
            # Create Series of PASSED variants, failed variants and non-variant positions
            # prioritising low coverage > PASSED variants > failed variants > non-variant (ref)
            all_vars = pd.concat([low_cov_vars, snpSeries, ufSeries], axis=1)
            # Remove duplicate Contig:Pos:Alt columns
            all_vars = all_vars.loc[:,~all_vars.columns.duplicated(keep='first')]

            # Rename columns with just the position
            # var = float(pos.alt), int(var) == pos
            new_colnames = [int(var) for var in all_vars.columns]
            all_vars.columns = new_colnames

            all_vars_dups = all_vars.loc[:,all_vars.columns.duplicated(keep=False)]
            all_vars = all_vars.loc[:,~all_vars.columns.duplicated(keep=False)]
            all_vars_dups = all_vars_dups.groupby(all_vars_dups.columns, axis=1).apply(
            lambda x: x.apply(
            lambda y: IUPAC_CODE[int(y[~y.isna()].astype('str').str.cat())], axis=1))

            all_vars_dups = all_vars_dups.astype('uint8')
            all_vars = all_vars.astype('uint8')

            all_vars = pd.concat([all_vars, all_vars_dups, ref_pos_series], axis=1)
            all_vars = all_vars.loc[:,~all_vars.columns.duplicated(keep='first')]

            return all_vars

        else: # no SNPs in contig, only coverage
            new_colnames = [int(var) for var in low_cov_vars.columns]
            low_cov_vars.columns = new_colnames
            low_cov_vars = low_cov_vars.loc[:,~low_cov_vars.columns.duplicated(keep='first')]

            return low_cov_vars

    def getLowCov(self, var_list, contig):
        """ Return series with each SNP position (column) in posList populated
        with low coverage """
        low_cov_cols = []
        pos_list = [int(pos) for pos in var_list]
        pos_idx = np.where(np.isin(pos_list,self.covSeries.get(contig)))
        low_cov_cols = [var_list[i] for i in pos_idx[0]]
        # low_cov_cols = [var_list[e] for e,x in enumerate(varList) if x.pos in covCols]

        return low_cov_cols
        # return self.covSeries.get(contig)

    def readPosFile(self, pos_file):
        with open(pos_file, "r") as pos_list:
            for line in pos_list:
                self.pos.append(int(line.strip()))

    def writeVcf(self):
        """ Writes out all variants in PASS and FAIL vcf """
        outfile = open("{}{}_freebayes_FILTERED.vcf".format(self.outdir,self.name), "w")
        outfile.write(self.header)
        for contig in self.variants.keys():
            for pos in self.variants[contig].keys():
                outfile.write(self.variants[contig][pos])
        outfile.close()

    def writeSampleVcf(self, posList, filename):
        """ Writes out vcf file if variants occurred based on list of given positions."""
        outfile = open(filename, "w")
        outfile.write(self.header)
        for contig in posList.keys():
            for pos in sorted(posList[contig]):
                var = self.variants[contig].get(pos)
                if var:
                    outfile.write(var)
        outfile.close()

# TODO: Load in the json format into dataframes, converting string to int
class jsonSample(Sample):
    """ Load Sample data from json file """
    def __init__(self,name,outdir,frac,filter_clust,filter_cliffs,contigs,header,refsUnfilt,snpSeries,ufSeries,covSeries,cliffs,mix,clusters):
        self.name = name
        self.outdir = outdir
        self.frac = frac
        self.filter_clust = filter_clust
        self.filter_cliffs = filter_cliffs
        self.contigs = contigs
        self.header = header
        self.refsUnfilt = refsUnfilt
        self.snpSeries = self.json_to_dataframe(snpSeries)
        self.ufSeries = self.json_to_dataframe(ufSeries)
        self.covSeries = covSeries
        self.cliffs = cliffs
        self.mix = mix
        self.clusters = clusters

    def json_to_dataframe(self, json_data):
        attribute = {}
        for key in json_data.keys():
            attribute[key] = pd.read_json(json_data[key]).astype('Int8')
        return attribute

# JSON Encoder for Sample class
class SampleEncoder(JSONEncoder):
    def default(self, object):
        if isinstance(object, Sample):
            return object.__dict__
        elif isinstance(object, set):
            return list(object)
        elif isinstance(object, pd.DataFrame):
            try:
                return object.to_json()
            except ValueError as e:
                print(object.columns[object.columns.duplicated(keep=False)])
                raise e
        else:
            # raises exception for unsupported types
            print("Not instance of Sample")
            return json.JSONEncoder.default(self, object)

def best_fit_slope_and_intercept(x,y):
    m = ((mean(x)*mean(y)) - mean(x*y)) / ((mean(x)*mean(x)) - mean(x*x))
    b = mean(y) - m*mean(x)
    return m, b

def squared_error(y_orig,y_line):
    return sum((y_line - y_orig) * (y_line - y_orig))

def coefficient_of_determination(y_orig,y_line):
    y_mean_line = [mean(y_orig) for y in y_orig]
    squared_error_regr = squared_error(y_orig, y_line)
    squared_error_y_mean = squared_error(y_orig, y_mean_line)
    return 1 - (squared_error_regr/squared_error_y_mean)

def read_cigar(cig_string):
    CIGAR = []
    prev = None
    l = None # length in cigar
    for c in cig_string:
        if prev is None:
            if c.isnumeric():
                l = c
            else:
                print(cig_string)
                print("Malformed CIGAR")
                exit()
        elif prev.isalpha():
            if c.isnumeric():
                l = c
            else:
                print(cig_string)
                print("Malformed CIGAR")
                exit()
        elif prev.isnumeric():
            if c.isnumeric(): # e.g. CIGAR = 13M
                l += c
            else:
                CIGAR.append([int(l),c])
        prev = c
    return CIGAR
