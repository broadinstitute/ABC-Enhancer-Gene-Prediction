import pandas as pd
import numpy as np
import os
import os.path
from tools import *
from neighborhoods import *

def make_candidate_regions_from_summits(macs_peaks, accessibility_file, genome_sizes, regions_includelist, regions_blocklist, n_enhancers, peak_extend, outdir):
    ## Generate enhancer regions from MACS summits: 1. Count reads in DHS peaks 2. Take top N regions, get summits, extend summits, merge
    
    outfile = os.path.join(outdir, os.path.basename(macs_peaks) + ".candidateRegions.bed") # output file name for candidateRegions.bed
    raw_counts_out = [] # initialize list for output file names
    for access_in in accessibility_file: # loop through input accessibilty files
       raw_counts_out.append(os.path.join(outdir, os.path.basename(macs_peaks) + "." + os.path.basename(access_in) + ".Counts.bed")) # add corresponding output file name
    if regions_includelist:
    	includelist_command = "(bedtools intersect -a {regions_includelist} -b {genome_sizes}.bed -wa | cut -f 1-3 && cat) |"
    else:
    	includelist_command = ""

    if regions_blocklist:
    	blocklist_command = "bedtools intersect -v -wa -a stdin -b {regions_blocklist} | "
    else:
    	blocklist_command = ""

    #1. Count DHS/ATAC reads in candidate regions
    #run_count_reads(accessibility_file, raw_counts_outfile, macs_peaks, genome_sizes, use_fast_count=True)
    
    # 1. Count DHS/ATAC reads in candidate regions for all accessibility files provided, and return the filename of the average # reads
    reads_out = count_reads_over_peaks(accessibility_file, raw_counts_out, macs_peaks, genome_sizes, outdir, use_fast_count=True)

    #2. Take top N regions, get summits, extend summits, merge, remove blocklist, add includelist, sort and merge
    #use -sorted in intersect command? Not worth it, both files are small
    command = "bedtools sort -i {reads_out} -faidx {genome_sizes} | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n {n_enhancers} |" + \
        "bedtools intersect -b stdin -a {macs_peaks} -wa |" + \
        "awk '{{print $1 \"\\t\" $2 + $10 \"\\t\" $2 + $10}}' |" + \
        "bedtools slop -i stdin -b {peak_extend} -g {genome_sizes} |" + \
        "bedtools sort -i stdin -faidx {genome_sizes} |" + \
        "bedtools merge -i stdin | " + \
        blocklist_command + \
        "cut -f 1-3 | " + includelist_command + \
        "bedtools sort -i stdin -faidx {genome_sizes} | bedtools merge -i stdin > {outfile}"

    command = command.format(**locals())

    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    print("Running: " + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')

    return stdoutdata


def make_candidate_regions_from_peaks(macs_peaks, accessibility_file, genome_sizes, regions_includelist, regions_blocklist, n_enhancers, peak_extend, minPeakWidth, outdir):
    ## Generate enhancer regions from MACS narrowPeak - do not use summits
    outfile = os.path.join(outdir, os.path.basename(macs_peaks) + ".candidateRegions.bed")
    raw_counts_out = [] # initialize list for output file names
    for access_in in accessibility_file: # loop through input accessibilty files
       raw_counts_out.append(os.path.join(outdir, os.path.basename(macs_peaks) + "." + os.path.basename(access_in) + ".Counts.bed")) # add corresponding output file name
       
    if regions_includelist:
        includelist_command = "(bedtools intersect -a {regions_includelist} -b {genome_sizes}.bed -wa | cut -f 1-3 && cat) |"
    else:
        includelist_command = ""

    if regions_blocklist:
        blocklist_command = "bedtools intersect -v -wa -a stdin -b {regions_blocklist} | "
    else:
        blocklist_command = ""

    #1. Count DHS/ATAC reads in candidate regions
    reads_out = count_reads_over_peaks(accessibility_file, raw_counts_out, macs_peaks, genome_sizes, outdir, use_fast_count=True)

    #2. Take top N regions, extend peaks (min size 500), merge, remove blocklist, add includelist, sort and merge
    #use -sorted in intersect command? Not worth it, both files are small
    command = "bedtools sort -i {reads_out} -faidx {genome_sizes} | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n {n_enhancers} |" + \
        "bedtools intersect -b stdin -a {macs_peaks} -wa |" + \
        "bedtools slop -i stdin -b {peak_extend} -g {genome_sizes} |" + \
        "awk '{{ l=$3-$2; if (l < {minPeakWidth}) {{ $2 = $2 - int(({minPeakWidth}-l)/2); $3 = $3 + int(({minPeakWidth}-l)/2) }} print $1 \"\\t\" $2 \"\\t\" $3}}' |" + \
        "bedtools sort -i stdin -faidx {genome_sizes} |" + \
        "bedtools merge -i stdin | " + \
        blocklist_command + \
        "cut -f 1-3 | " + includelist_command + \
        "bedtools sort -i stdin -faidx {genome_sizes} | bedtools merge -i stdin > {outfile}"

    command = command.format(**locals())

    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    print("Running: " + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')
    if not err == '':
        raise RuntimeError("Command failed.")

    return stdoutdata
    
# count reads over however many DHS files and return average
def count_reads_over_peaks(accessibility_file, raw_counts_out, macs_peaks, genome_sizes, outdir, use_fast_count=True):
    for access_in, counts_out in zip(accessibility_file, raw_counts_out):
        run_count_reads(access_in, counts_out, macs_peaks, genome_sizes, use_fast_count)
        nFiles = len(raw_counts_out)
        
    if nFiles > 1:
        avg_out = os.path.join(outdir, os.path.basename(macs_peaks) + ".averageAccessibility.Counts.bed")
        df1 = pd.read_csv(raw_counts_out[0], sep='\t', header=None)
        totalCounts = df1[3]
        for i in range(1, nFiles):
            dfx = pd.read_csv(raw_counts_out[i], sep='\t', header=None, usecols=[3])
            print(len(dfx))
            df1[3] = df1[3].add(dfx[3])
        df1[3] = df1[3]/nFiles
        df1.to_csv(avg_out, header=None, index=None, sep="\t")
        return avg_out
    else:
        return raw_counts_out[0]
