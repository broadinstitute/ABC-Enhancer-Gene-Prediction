import os
import numpy as np
import pandas as pd
import re
from subprocess import check_call
import sys
import pyranges as pr

# setting this to raise makes sure that any dangerous assignments to pandas
# dataframe slices/subsets error rather than warn
pd.set_option('mode.chained_assignment', 'raise')

def run_command(command, **args):
    print("Running command: " + command)
    return check_call(command, shell=True, **args)

def write_connections_bedpe_format(pred, outfile, score_column):
    #Output a 2d annotation file with EP connections in bedpe format for loading into IGV
    pred = pred.drop_duplicates()

    towrite = pd.DataFrame()

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

def determine_expressed_genes(genes, expression_cutoff, activity_quantile_cutoff):
    #Evaluate whether a gene should be considered 'expressed' so that it runs through the model

    #A gene is runnable if:
    #It is expressed OR (there is no expression AND its promoter has high activity)

    genes['isExpressed'] = np.logical_or(genes.Expression >= expression_cutoff, np.logical_and(np.isnan(genes.Expression), genes.PromoterActivityQuantile >= activity_quantile_cutoff))

    return(genes)

def write_params(args, file):
    with open(file, 'w') as outfile:
        for arg in vars(args):
            outfile.write(arg + " " + str(getattr(args, arg)) + "\n")

def df_to_pyranges(df, start_col='start', end_col='end', chr_col='chr', start_slop=0, end_slop=0):
    df['Chromosome'] = df[chr_col]
    df['Start'] = df[start_col] - start_slop
    df['End'] = df[end_col] + end_slop

    return(pr.PyRanges(df))
