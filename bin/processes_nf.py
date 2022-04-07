#!/usr/bin/env python3


import csv
import numpy as np
import pandas as pd
import pyranges as pyr
from pysam import VariantFile
import scipy.sparse as sp
from statistics import mean
import multiprocessing
from multiprocessing.pool import ThreadPool
import os
import sys
import datetime
from Bio import SeqIO
# import tracemalloc
from dna_dicts import *
from sample import *

# Globals
MIN_DEPTH = 10
SNP_CLUST = 3
WINDOW = 10
ALT_RATIO = 0.75
CLIFF_WINDOW = 10
CLIFF_SLOPE = 3 # 1 for testing, change to 3 for running
SEARCH_WINDOW = 10
R_SQUARED_THRESHOLD = 0.7
MIN_MQM = 30

class Results:

    def __init__(self,  min_cov, samples = [], excluded = []):
        self.min_cov = min_cov
        self.samples = samples # List of class Sample (initialise with list of Sample or addSample())
        self.excluded = excluded # samples excluded due to low coverage %

    def addSample(self, sample):
        if self.checkExclude(sample):
            self.excluded.append(sample)
        else:
            self.samples.append(sample)

    def addNewSample(self, sample_name, frac, cov_file, vcf_file):
        # add sample from vcf and coverage file
        sample = Sample(sample_name, frac, cov_file, vcf_file)
        if self.checkExclude(sample):
            self.excluded.append(sample)
        else:
            self.samples.append(sample)

    def checkExclude(self, sample):
        exclude = False
        for contig in sample.coverage.keys():
            if sample.coverage[contig] < self.min_cov:
                exclude = True
        return exclude

    def getAllFiltered(self):
        """ Returns dict{contig:set{pos.alt}} of all high quality variants and low coverage positions across all samples 
        Where pos.base is a float """
        allVars = {}
        for sample in self.samples:
            vars = sample.getFiltered()
            for contig in vars.keys():
                allVars.setdefault(contig,set()).update(vars[contig].columns)
        return allVars

    def getAllUnfiltered(self, contig, ref_series, ref_pos_series):
        """ Returns list of pandas sub-series of unfiltered variants and low coverage for each sample based on colnames for a particular reference contig """
        vars = []
        for sample in self.samples:
            vars.append(sample.getUnfiltered(contig, ref_series, ref_pos_series))
        return vars

    def getSampleUnfiltered(self, args):
        sample, contig, ref_series, ref_pos_series = args
        return sample.getUnfiltered(contig, ref_series, ref_pos_series)

    def getAllLowCov(self, var_list, contig):
        low_cov = set()
        for sample in self.samples:
            low_cov.update(sample.getLowCov(var_list,contig))
        return low_cov

    def getSnpRef(self):
        """ Returns dictionary of snp positions and reference alleles (key=CONTIG:POS) """
        ref = {}
        for sample in self.samples:
            for contig in sample.refsUnfilt.keys():
                if ref.get(contig) is None:
                    ref[contig] = sample.refsUnfilt[contig] # Dict keys = str("Contig:Pos")
                else:
                    ref[contig].update(sample.refsUnfilt[contig])
        return ref

    def writeAllVcf(self, outdir, posList=None, suffix="filtered"):
        """ Writes out vcf file if variants occurred based on list of given positions.
        If no pos list given, all variants are written out."""
        for sample in self.samples:
            plist = posList
            filename = outdir+sample.name+"_"+suffix+".vcf"
            if plist is None:
                plist = {}
                for contig in sample.variants.keys():
                    plist[contig] = sample.variants.get(contig).keys()
            sample.writeSampleVcf(plist, filename)

    def getCliffs(self):
        cliffs = {}
        for sample in self.samples:
            for contig in sample.contigs:
                if cliffs.get(contig) is None:
                    cliffs[contig] = set()
                if sample.cliffs.get(contig):
                    cliffs[contig] = cliffs[contig].union(sample.cliffs.get(contig))
        return cliffs

    def getMix(self):
        mix = {}
        for sample in self.samples:
            for contig in sample.contigs:
                if mix.get(contig) is None:
                    mix[contig] = set()
                if sample.mix.get(contig):
                    mix[contig] = mix[contig].union(sample.mix.get(contig))
        return mix

    def getClusters(self):
        clusters = {}
        for sample in self.samples:
            for contig in sample.contigs:
                if clusters.get(contig) is None:
                    clusters[contig] = set()
                if sample.clusters.get(contig):
                    clusters[contig] = clusters[contig].union(sample.clusters.get(contig))
        return clusters

class Mtx:
    """ SNP and coverage matrices
    Attr:
    results (Results)
    self.mtx (Pandas DataFrame): All SNP call results and low coverage regions
    """

    def __init__(self, ref, sample, results, maskFile = None, threads = 1, recomb = True, exclude = False, outdir=""):
        self.ref = ref
        self.sample = sample
        self.contigs = set()
        self.ref_lists = {}
        self.ref_dicts = {}
        self.results = results
        self.threads = threads
        self.mask_file = maskFile
        self.recomb = recomb
        self.exclude = exclude
        self.outdir = outdir
        self.mtx = {} # all snps and all low coverage positions
        self.snp = {} # all snps
        self.core = {} # all core snps
        self.non_core_pos = {} # {contig:[int,int]}

        
        self.__createRefFrame()

        if self.sample not in results.excluded:
            self.__buildMtx() # Matrix of all variant and low coverage positions
            self.__getSnpMtx()
            self.__getCoreMtx()

            print("Writing full fastas files:", datetime.datetime.now())
            self.writeFastaAlns()

            print("Writing fuzzy SNP alignments:", datetime.datetime.now())
            self.writeFullSnps()

            print("Writing core SNP alignments:", datetime.datetime.now())
            self.writeCoreSnps()

        else:
            self.writeExcludedStats()



    def __createRefFrame(self):
        """ Load in reference as dataframe """
        for seq_record in SeqIO.parse(self.ref, 'fasta'):
            self.contigs.add(seq_record.id)
            vals = list('-'+seq_record.seq) # add leading N to turn indexing into 1-based
            self.ref_lists[seq_record.id] = vals
            # create dictionary of reference alleles for every position
            self.ref_dicts[seq_record.id] = { i : vals[i] for i in range(0, len(vals))}

    def __buildMtx(self):
        """ Build pandas matrix of low coverage regions and all snps where at least one snp call was 'high quality' """
        # Get all low coverage positions and high quality variants (columns) for alignment/matrix
        filtVars = self.results.getAllFiltered()

        for contig in filtVars.keys():
            ref_dict_sub = dict(map(lambda key: (key, DNA_CODE[self.ref_dicts[contig].get(int(key))]), sorted(list(filtVars[contig]))))
            # Create reference series for every filtered position
            # colnames of series are class sample.SNP
            ref_series = pd.DataFrame(ref_dict_sub, index=[0], dtype='uint8')

            # Prepare ref series for Concatenating
            # Need to remove duplicate variant positions e.g. contig:1:A, contig:1:C
            # by renaming columns to the position and removing duplicates
            ref_pos_series = ref_series.copy(deep=True)
            ref_colnames = [int(var) for var in ref_series.columns]
            ref_pos_series.columns = ref_colnames
            ref_pos_series = ref_pos_series.loc[:,~ref_pos_series.columns.duplicated(keep='first')]

            # Single sample version
            # Returns series with colnames as position
            snps = self.results.getSampleUnfiltered([self.sample, contig, ref_series.copy(deep=False), ref_pos_series.copy(deep=False)])

            if snps is not None: # snps in contig
                self.snp[contig] = snps

    def merge_cols(self, y):
        if y.isna().all():
            return pd.NA
        else:
            return(IUPAC_CODE[int(y[~y.isna()].astype('str').str.cat())])

    def __getSnpMtx(self):
        """
        Return full SNP matrix removing recombination, mixed alleles and those that occur in cliffs
        Removes sites describing coverage only (i.e. keep only columns where at least 1 values is A,T,G,C)
        """
        # refs = self.results.getSnpRef()

        cliffs = self.results.getCliffs()
        # self.results.writeAllVcf(self.outdir, cliffs, "cliffs")

        mix = self.results.getMix()
        # self.results.writeAllVcf(self.outdir, mix, "mix")

        if self.recomb:
            clusters = self.results.getClusters()
            # self.results.writeAllVcf(self.outdir, clusters, "clusters")
        else:
            clusters = None # empty set if not using inbuilt rough recombination removal
        if self.mask_file:
            masked = self.maskedPos()
        else:
            masked = None

        for contig in self.snp.keys():
            cliff_pos = cliffs.get(contig)
            if mix is not None:
                mix_pos = mix.get(contig)
            else:
                mix_pos = []
            if clusters is not None:
                clust_pos = clusters.get(contig)
            else:
                clust_pos = []
            if masked is not None:
                mask_pos = masked.get(contig)
            else:
                mask_pos = []

            print("Number of cliff positions for {}: {}".format(contig, len(cliff_pos)))
            print("Number of mix/ambigious SNP call positions for {}: {}".format(contig, len(mix_pos)))
            print("Number of clustered positions for {}: {}".format(contig, len(clust_pos)))
            print("Number of masked positions from bed/gff file {}: {}".format(contig, len(mask_pos)))

            # Set as argument to remove all ambiguous positions, default = keep
            # remove = set().union(cliff_pos, mask_pos, clust_pos, mix_pos)
            remove = set().union(cliff_pos, mask_pos, clust_pos)

            print("Total number of all masked/removed positions (cliffs, recombination and bed/gff file) for {}: {}".format(contig, len(remove)))

            keep = list(set(self.snp[contig].columns).difference(remove))
            self.snp[contig] = self.snp[contig][keep]

            # # Removes columns where all values are the same
            # # This will exclude SNPs were all samples have them relative to the reference
            # # only if -e (exclude) is specified to exclude reference from alignment
            # # SELECT TO INCLUDE OR EXCLUDE IN PROD VERSION
            # nuniq = self.snp[contig].nunique()
            # to_drop = nuniq[nuniq==1].index
            # print("Dropping non-unique variants",datetime.datetime.now())
            # self.snp[contig] = self.snp[contig].drop(to_drop, axis=1)

            print("Reindexing",datetime.datetime.now())
            self.snp[contig] = self.snp[contig].reindex(sorted(self.snp[contig].columns), axis=1)

            self.snp.get(contig).replace(IUPAC_STR).to_csv("{}_{}_full_snps.csv".format(self.sample.name, contig), index_label="SAMPLE")
        
            print("Done",datetime.datetime.now())

    def __getCoreMtx(self):
        """ Return core SNP matrix (i.e. remove all positions with 0 low coverage in any sample) """
        # core_pos = {}
        for contig in self.snp.keys():
            # Positions where all samples have coverage
            keep = sorted(list(set(self.snp[contig].columns).difference(self.results.getAllLowCov(self.snp[contig].columns, contig))))

            self.core[contig] = self.snp[contig][keep]

            self.core.get(contig).replace(IUPAC_STR).to_csv("{}_{}_core_snps.csv".format(self.sample.name, contig), index_label="SAMPLE")
            

    def getPairwiseDist(self):
        """ Calculate pairwise counts of core SNP differences """
        pw = None
        for contig in self.core.keys():
            if pw is None:
                pw = np.zeros((len(self.core[contig].index),len(self.core[contig].index)))
            pw += ((self.core[contig].values ^ self.core[contig].values[:,None]) > 0).sum(2)
        pw_mtx = pd.DataFrame(pw, index=self.core[contig].index.to_list(), columns=self.core[contig].index.to_list(), dtype='uint16')

        return pw_mtx

    def maskedPos(self):
        """ Read in bed and/or gff file(s) of positions to mask.
        Bed format: 0-based coordinates
        chrom   start   end
        GFF format: 1-based coordinates
        chrom   source  feature start   end score   strand  frame   group/info
        Returns:
            """
        try:
            ranges = pyr.read_gff3(self.mask_file)
            file_format = "gff3"
        except:
            try:
                ranges = pyr.read_gff(self.mask_file)
                file_format = "gff"
            except:
                try:
                    ranges = pyr.read_bed(self.mask_file)
                    file_format = "bed"
                except Exception as e:
                    print(e)
                    sys.exit("Unable to parse gff/bed file.")

        print("Loaded {} regions for masking from {} formated file".format(len(ranges), file_format))
        print("PYRANGES CHROMOSOMES", ranges.chromosomes)
        posRanges = pyr.PyRanges(chromosomes=ranges.chromosomes[0], starts=pos, ends=pos)
        maskRanges = posRanges.intersect(ranges)
        return(set(maskRanges.Start))

    def concatCoreMtx(self):
        """ concat all core matrices for each contig """
        concat_mtx = pd.DataFrame(index = [self.sample.name])
        for contig in self.core.keys():
            mtx = self.core[contig]
            # mtx = mtx.rename(index={contig:'reference'})
            concat_mtx = pd.concat([concat_mtx, mtx], axis=1, sort=False) # Full outer join
        return(concat_mtx)

    def concatFullMtx(self):
        """ concat all full snp matrices for each contig """
        concat_mtx = pd.DataFrame(index = [self.sample.name])
        for contig in self.snp.keys():
            mtx = self.snp[contig]
            # Replace NA with reference allele
            mtx = mtx.fillna(self.ref_dicts[contig][0]) # inplace=True very slow
            # mtx = mtx.rename(index={contig:'reference'})
            concat_mtx = pd.concat([concat_mtx, mtx], axis=1, sort=False)

        return(concat_mtx)

    def writeFastaAlns(self):
        """
        Generate fasta files for each sample
        """
        print("Writing Fastas:", datetime.datetime.now())
        stats_out = open("{}_core_stats.tsv".format(self.sample.name), 'w')
        snps_out = open("{}_snpmatrix.csv".format(self.sample.name), 'w')
        snps_out.write("Contig,POS,{}\n".format(self.sample.name))
        full_out = open("{}_full_aln.fasta".format(self.sample.name), 'w')

        full_out.write(">{}\n".format(self.sample.name))

        for contig in self.snp.keys():
            full_length = len(self.ref_lists[contig]) -1 # -1 due to leading additional '-'

            # single_sample
            self.writeFastaHelper([self.sample.name, contig, full_length, full_out, stats_out, snps_out])

        full_out.write("\n")
        full_out.close()

        stats_out.close()

    def writeFastaHelper(self, args):
        sample, contig, full_length, full_out, stats_out, snps_out = args

        print("Writing fasta for {} {}".format(contig,sample), datetime.datetime.now())

        print("Generate sample_series for {} {}".format(contig,sample), datetime.datetime.now())
        sample_series = self.snp[contig].loc[sample]

        print("Replace with IUPAC_STR for {} {}".format(contig,sample), datetime.datetime.now())
        sample_series = sample_series.replace(IUPAC_STR, inplace=False)

        print("Getting reference list".format(contig,sample), datetime.datetime.now())
        fasta_list = self.ref_lists[contig].copy()

        print("Replacing values".format(contig,sample), datetime.datetime.now())
        for (index, replacement) in zip(sample_series.index, sample_series.values):
            fasta_list[index] = replacement

        print("Writing out snps for {} {}".format(contig,sample), datetime.datetime.now())
        sample_df = sample_series.to_frame()
        sample_df.insert(0, "POS", sample_series.index, True)
        sample_df.insert(0, "Contig", [contig] * len(sample_series.index), True)
        sample_df.to_csv(snps_out, sep=',', index=False, header=False)

        print("Writing to file".format(contig,sample), datetime.datetime.now())
        full_out.write("".join(fasta_list[1:]))

        if self.sample.covSeries.get(contig) is not None:
            core_length = full_length - len(self.sample.covSeries[contig])
        else:
            core_length = full_length
        
        stats_string = "{}\t{}\t{}\t{}\n".format(sample, contig, core_length, round(core_length/full_length*100,2))
        stats_out.write(stats_string)

    def writeExcludedStats(self):
        stats_out = open("{}_excluded.tsv".format(self.sample.name), 'w')
        for contig in self.sample.coverage.keys():
            out_string = "{}\t{}\t{}\n".format(self.sample.name, contig, self.sample.coverage[contig])
            stats_out.write(out_string)
        stats_out.close()
        

    # TODO: Speed up using lambda https://towardsdatascience.com/how-to-make-your-pandas-loop-71-803-times-faster-805030df4f06
    def writeFullSnps(self):
        """ Write out fasta alignment files """
        # Whole genome SNP alignment
        mtx = self.concatFullMtx()
        mtx = mtx.replace(IUPAC_STR, inplace=False)
        with open("{}_full_snp.fasta".format(self.sample.name), 'w') as out:
            for e, sample in enumerate(mtx.index):
                val_list = list(mtx.iloc[e, :])
                out.write(">{}\n".format(sample))
                try:
                    out.write("".join(val_list))
                except:
                    print("Error writing full SNP Fasta",val_list)
                out.write("\n")
        out.close()

    def writeCoreSnps(self):
        """ Write out core fasta and phylip """
        mtx = self.concatCoreMtx()
        mtx = mtx.replace(IUPAC_STR, inplace=False)
        with open("{}_core_snp.fasta".format(self.sample.name), 'w') as out:
            for e, sample in enumerate(mtx.index):
                val_list = list(mtx.iloc[e, :])
                out.write(">{}\n".format(sample))
                out.write("".join(val_list))
                out.write("\n")
        out.close()

        # # Phylip file
        # num_samples = len(mtx.index)
        # total_snps = len(mtx.columns)
        # mtx.insert(0,0,np.nan) # whitespace after sample names
        # mtx.index = [name[0:10] for name in mtx.index] # phylip sample name convention
        #
        # # start = llength+1
        # # # break down into blocks/alignments of length llength
        # # phy += mtx.iloc[:,0:start].to_csv(na_rep=' ', sep=',', header=False).replace(',', '')
        # # print(phy)
        # # for n in range(start+llength,num_snps,llength):
        # #     phy+=mtx.iloc[:,start:n].to_csv(na_rep=' ', sep=',', index=None, header=None).replace(',', '')
        # #     start = n
        # # if start <= num_snps:
        # #     phy+=mtx.iloc[:,start:].to_csv(na_rep=' ', sep=',', index=None, header=None).replace(',', '')
        #
        # phy_filename = "{}core.phy".format(self.outdir)
        # phy_file = open(phy_filename, "w")
        # header = str(num_samples)+" "+str(total_snps)+"\n"
        # phy_file.write(header)
        # phy = mtx.to_csv(na_rep=' ', sep=',', header=None).replace(',','')
        # phy_file.write(phy)
        # phy_file.close()

# Additional functions
def remove_from_list_helper(args):
    sample, pos_list, core_fasta = args
    all_pos = [n for n in range(0,len(core_fasta))]
    all_pos = set(all_pos)
    keep_pos = sorted(list(all_pos.difference(set(pos_list))))
    keep_fasta = [core_fasta[p] for p in keep_pos]
    return (sample,keep_fasta)

def load_results(samples, outdir):
    """ Read in bcf and vcf files for each sample into Results """
    results = Results()
    for sample in samples:
        cov_file = outdir+sample+'.pileup'
        vcf_file = outdir+sample+'.bcf'
        results.addSample(sample, frac, cov_file, vcf_file)
    return results

def load_sample(args_list):
    sample_name, frac, cov_file, vcf_file = args_list
    print("Loading",sample_name, "on", os.getpid())
    sample = Sample(sample_name, frac, cov_file, vcf_file)
    return sample

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
