#!/usr/bin/env python3

import os
import pandas
import pickle
from intervaltree import IntervalTree, Interval
import pysam
import numpy as np
import pandas as pd
import re
from subprocess import check_call
import sys

# setting this to raise makes sure that any dangerous assignments to pandas
# dataframe slices/subsets error rather than warn
pandas.set_option('mode.chained_assignment', 'raise')


def read_genes(filename):
    genes = pandas.read_table(filename)

    #Deduplicate genes.
    gene_cols = ['chr','tss','name','Expression','PromoterActivityQuantile']
    genes = genes[gene_cols]
    genes.drop_duplicates(inplace=True)

    return genes


def get_gene_name(gene):
    try:
        out_name = gene['name'] #if ('symbol' not in gene.keys() or gene.isnull().symbol) else gene['symbol']
    except:
        out_name = "_UNK"
    return str(out_name)

def get_score_filename(gene, outdir = None):
    out_name = get_gene_name(gene)
    outfile = "{}_{}_{}.prediction.txt.gz".format(out_name, gene.chr, int(gene.tss))
    if outdir is not None:
        outfile = os.path.join(outdir, outfile)
    return outfile


def write_scores(outdir, gene, enhancers):
    outfile = get_score_filename(gene)
    enhancers.to_csv(os.path.join(outdir, outfile), sep="\t", index=False, compression="gzip", float_format="%.6f", na_rep="NaN")


# class GenomicRangesIntervalTree(object):
#     def __init__(self, filename, slop=0, isBed=False):
#         if isBed:
#             self.ranges = pandas.read_table(filename, header=None, names=['chr', 'start', 'end', 'Score'])
#         else:
#             self.ranges = pandas.read_table(filename)

#         self.ranges['start'] = self.ranges['start'] - slop
#         self.ranges['end'] = self.ranges['end'] + slop
#         self.ranges['end'] = [ max(x,y) for x,y in zip(self.ranges['start']+1,self.ranges['end']) ]
#         assert(pandas.DataFrame.all(self.ranges.start <= self.ranges.end))
        
#         self.intervals = {}
#         for chr, chrdata in self.ranges.groupby('chr'):
#             self.intervals[chr] = IntervalTree.from_tuples(zip(chrdata.start,
#                                                                chrdata.end,
#                                                                chrdata.index))

#     def within_range(self, chr, start, end):
#         # Returns empty data frame (0 rows) if there is no overlap
#         if start == end:   ## Interval search doesn't like having start and end equal
#             end = end + 1
#         result = self.ranges.iloc[[], :].copy()
#         if chr in self.intervals:
#             overlaps = self.intervals[chr][start:end]
#             indices = [idx for l, h, idx in overlaps]
#             result = self.ranges.iloc[indices, :].copy()
#         return result

#     def overlaps(self, chr, start, end):
#         locs = self.intervals[chr].overlaps()

#     def __getitem__(self, idx):
#         return self.ranges[idx]

class GenomicRangesIntervalTree(object):
    def __init__(self, filename, slop=0, isBed=False):
        if isBed:
            self.ranges = pandas.read_table(filename, header=None, names=['chr', 'start', 'end', 'Score'])
        elif isinstance(filename, pandas.DataFrame):
            self.ranges = filename
        else:
            self.ranges = pandas.read_table(filename)

        self.ranges['start'] = self.ranges['start'] - slop
        self.ranges['end'] = self.ranges['end'] + slop
        self.ranges['end'] = [ max(x,y) for x,y in zip(self.ranges['start']+1,self.ranges['end']) ]
        assert(pandas.DataFrame.all(self.ranges.start <= self.ranges.end))
        
        self.intervals = {}
        for chr, chrdata in self.ranges.groupby('chr'):
            self.intervals[chr] = IntervalTree.from_tuples(zip(chrdata.start,
                                                               chrdata.end,
                                                               chrdata.index))

    def within_range(self, chr, start, end):
        # Returns empty data frame (0 rows) if there is no overlap
        if start == end:   ## Interval search doesn't like having start and end equal
            end = end + 1
        result = self.ranges.iloc[[], :].copy()
        if chr in self.intervals:
            overlaps = self.intervals[chr][start:end]
            indices = [idx for l, h, idx in overlaps]
            result = self.ranges.iloc[indices, :].copy()
        return result

    def overlaps(self, chr, start, end):
        locs = self.intervals[chr].overlaps()

    def __getitem__(self, idx):
        return self.ranges[idx]

def read_enhancers(filename):
    return GenomicRangesIntervalTree(filename)

class DataCache(object):
    def __init__(self, directory):
        os.makedirs(directory, exist_ok=True)
        self.directory = directory

    def __contains__(self, filename):
        cache_name = os.path.join(self.directory, filename.replace(os.sep, '__'))
        if os.path.exists(cache_name) and (os.path.getctime(cache_name) > os.path.getctime(filename)):
            return True
        return False

    def __getitem__(self, filename):
        cache_name = os.path.join(self.directory, filename.replace(os.sep, '__'))
        if os.path.exists(cache_name) and (os.path.getctime(cache_name) > os.path.getctime(filename)):
            with open(cache_name, "rb") as f:
                return pickle.load(f)
        raise KeyError

    def __setitem__(self, filename, value):
        cache_name = os.path.join(self.directory, filename.replace(os.sep, '__'))
        with open(cache_name, "wb") as f:
            pickle.dump(value, f, protocol=pickle.HIGHEST_PROTOCOL)

def run_command(command, **args):
    print("Running command: " + command)
    return check_call(command, shell=True, **args)

def write_connections_bedpe_format(pred, outfile, score_column):
    #Output a 2d annotation file with EP connections in bedpe format for loading into IGV
    pred = pred.drop_duplicates()

    towrite = pandas.DataFrame()

    towrite["chr1"] = pred["chr"]
    towrite["x1"] = pred['start']
    towrite["x2"] = pred['end']
    towrite["chr2"] = pred["chr"]
    towrite["y1"] = pred["TargetGeneTSS"]
    towrite["y2"] = pred["TargetGeneTSS"]
    towrite["name"] = pred["TargetGene"] + "_" + pred["name"]
    towrite["score"] = pred[score_column]
    towrite["strand1"] = "."
    towrite["strand2"] = "."

    towrite.to_csv(outfile, header=False, index=False, sep = "\t")

def check_gene_for_runnability(gene, expression_cutoff, activity_quantile_cutoff):
    #Evaluate whether a gene should be considered 'expressed' so that it runs through the model

    #A gene is runnable if:
    #It is expressed OR (there is no expression AND its promoter has high activity)

    try:
        gene['expressed'] = gene['Expression'] >= expression_cutoff
    except:
        gene['expressed'] = np.NaN

    is_active = gene["PromoterActivityQuantile"] >= activity_quantile_cutoff
    missing_expression = np.isnan(gene.Expression) or (gene.Expression is None) or (gene.Expression is "")
    should_run = (gene.expressed is True) or (missing_expression and is_active)

    return(should_run)
