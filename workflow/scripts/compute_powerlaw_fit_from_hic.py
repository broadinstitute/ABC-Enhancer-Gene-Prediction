import numpy as np
import sys, traceback
import pandas
from scipy import stats
import argparse
import glob
import os
from hic import *

#To do: 
#1. Use MLE to estimate exponent?
#2. Support bedpe

def parseargs():
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Helper to compute hic power-law fit parameters',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--hicDir', help="Directory containing observed HiC KR normalized matrices. File naming and structure should be: hicDir/chr*/chr*.KRobserved")
    parser.add_argument('--outDir', help="Output directory")
    parser.add_argument('--hic_type', default = 'juicebox', choices=['juicebox','bedpe'], help="format of hic files")
    parser.add_argument('--resolution', default=5000, type=int, help="For Juicebox: resolution of hic dataset (in bp). For bedpe: distances will be binned to this resolution for powerlaw fit")
    parser.add_argument('--minWindow', default=5000, type=int, help="Minimum distance between bins to include in powerlaw fit (bp). Recommended to be at least >= resolution to avoid using the diagonal of the HiC Matrix")
    parser.add_argument('--maxWindow', default=1000000, type=int, help="Maximum distance between bins to include in powerlaw fit (bp)")
    parser.add_argument('--chr', default='all', help="Comma delimited list of chromosomes to use for fit. Defualts to chr[1..22],chrX")

    args = parser.parse_args()
    return(args)

def main():
    args = parseargs()
    os.makedirs(args.outDir, exist_ok=True)

    HiC = load_hic_for_powerlaw(args)

    #Run 
    slope, intercept, hic_mean_var = do_powerlaw_fit(HiC)

    #print
    res = pandas.DataFrame({'resolution' : [args.resolution], 'maxWindow' : [args.maxWindow], 'minWindow' : [args.minWindow] ,'pl_gamma' : [slope], 'pl_scale' : [intercept] })
    res.to_csv(os.path.join(args.outDir, 'hic.powerlaw.txt'), sep='\t', index=False, header=True)

    hic_mean_var.to_csv(os.path.join(args.outDir, 'hic.mean_var.txt'), sep='\t', index=True, header=True)

def load_hic_for_powerlaw(args):
    if args.chr == 'all':
        chromosomes = ['chr' + str(x) for x in  list(range(1,23))] + ['chrX']
    else:
        chromosomes = args.chr.split(',')

    all_data_list = []
    for chrom in chromosomes:
        try:
            if args.hic_type == 'juicebox':
                hic_file, hic_norm_file, hic_is_vc = get_hic_file(chrom, args.hicDir, allow_vc=False)
                print("Working on {}".format(hic_file))
                this_data = load_hic(hic_file = hic_file, 
                    hic_norm_file = hic_norm_file,
                    hic_is_vc = hic_is_vc, 
                    hic_type = 'juicebox', 
                    hic_resolution = args.resolution, 
                    tss_hic_contribution = 100, 
                    window = args.maxWindow, 
                    min_window = args.minWindow, 
                    gamma = np.nan, 
                    interpolate_nan=False)
                this_data['dist_for_fit'] = abs(this_data['bin1'] - this_data['bin2']) * args.resolution
                all_data_list.append(this_data)
            elif args.hic_type == 'bedpe':
                hic_file, hic_norm_file, hic_is_vc = get_hic_file(chrom, args.hicDir, hic_type='bedpe')
                print("Working on {}".format(hic_file))
                this_data = load_hic(hic_file = hic_file,
                                     hic_type = 'bedpe',
                                     hic_norm_file = None,
                                     hic_is_vc = None, 
                                     hic_resolution = None, 
                                     tss_hic_contribution = None, 
                                     window = None, 
                                     min_window = None, 
                                     gamma = None)

                #Compute distance in bins as with juicebox data. 
                #This is needed to in order to maintain consistency, but is probably slightly less accurate.
                #Binning also reduces noise level.
                rawdist = abs((this_data['x2'] + this_data['x1'])/2 - (this_data['y2'] + this_data['y1'])/2)
                this_data['dist_for_fit'] = (rawdist // args.resolution) * args.resolution
                this_data = this_data.loc[np.logical_and(this_data['dist_for_fit'] >= args.minWindow, this_data['dist_for_fit'] <= args.maxWindow)]
                all_data_list.append(this_data)
            else:
                error('invalid --hic_type')

        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)

    all_data = pd.concat(all_data_list)

    return(all_data)

def do_powerlaw_fit(HiC):
    print("Running regression")

    #TO DO:
    #Print out mean/var plot of powerlaw relationship
    HiC_summary = HiC.groupby('dist_for_fit').agg({'hic_contact' : 'sum'})
    HiC_summary['hic_contact'] = HiC_summary.hic_contact / HiC_summary.hic_contact.sum() #technically this normalization should be over the entire genome (not just to maxWindow). Will only affect intercept though..
    res = stats.linregress(np.log(HiC_summary.index), np.log(HiC_summary['hic_contact']))

    hic_mean_var = HiC.groupby('dist_for_fit').agg({'hic_contact' : ['mean','var']})
    hic_mean_var.columns = ['mean', 'var']

    return res.slope, res.intercept, hic_mean_var

if __name__ == '__main__':
    main()
