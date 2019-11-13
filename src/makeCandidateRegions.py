import pandas as pd
import numpy as np
import argparse
import os
from peaks import *
import traceback
from tools import write_params

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Make peaks file for a given cell type',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    
    parser.add_argument('--narrowPeak', required=required_args, help="narrowPeak file output by macs2. Must include summits (--call-summits)")
    parser.add_argument('--bam', required=required_args, help="DNAase-Seq or ATAC-Seq bam file")
    parser.add_argument('--chrom_sizes', required=required_args, help="File listing chromosome size annotaions")
    parser.add_argument('--outdir', required=required_args)
    
    parser.add_argument('--nStrongestPeaks', default=175000, help="Number of peaks to use for defining candidate regions")
    parser.add_argument('--peakExtendFromSummit', default=250, help="Number of base pairs to extend each preak from its summit")

    parser.add_argument('--regions_whitelist', default="", help="Bed file of regions to forcibly include in candidate enhancers. Overrides regions_blacklist")
    parser.add_argument('--regions_blacklist', default="", help="Bed file of regions to forcibly exclude from candidate enhancers")
    
    args = parser.parse_args()
    return(args)

def processCellType(args):
	os.makedirs(os.path.join(args.outDir), exist_ok=True)
	write_params(args, os.path.join(args.outDir, "params.txt"))

	#Make candidate regions
	make_candidate_regions_from_summits(macs_peaks = args.narrowPeak, 
										accessibility_file = args.bam, 
										genome_sizes = args.chrom_sizes, 
										regions_whitelist = args.regions_whitelist,
										regions_blacklist = args.regions_blacklist,
										n_enhancers = args.nStrongestPeaks, 
										peak_extend = args.peakExtendFromSummit, 
										outdir = args.outDir)

def main(args):
    processCellType(args)

if __name__ == '__main__':
    args = parseargs()
    main(args)




