from hic import *
import numpy as np
from functools import reduce
import argparse
import sys, os, os.path
from tools import write_params
import pyranges

# To do
# Final output matrix needs to be KR normed as well
def parseargs():
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Make average HiC dataset',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--celltypes', required=True, help="Comma delimitted list of cell types")
    parser.add_argument('--chromosome', required=True, help="Chromosome to compute on")
    parser.add_argument('--basedir', required=True, help="Basedir")
    parser.add_argument('--outDir', required=True, help="Output directory")
    parser.add_argument('--resolution', default=5000, type=int, help="Resolution of hic dataset (in bp)")
    parser.add_argument('--ref_scale', default=5.41, type=float, help="Reference scale parameter")
    parser.add_argument('--ref_gamma', default=-.876, type=float, help="Reference gamma parameter")
    parser.add_argument('--min_cell_types_required', default=3, type=int, help="Minimum number of non-nan entries required to calculate average for a hic bin")

    args = parser.parse_args()
    return(args)

def main():
    args = parseargs()
    os.makedirs(args.outDir, exist_ok=True)

    #Write params file
    write_params(args, os.path.join(args.outDir, "params.txt"))

    #Parse cell types
    cell_types = args.celltypes.split(",")

    # chromosomes = ['chr' + str(x) for x in range(1,23)] + ['chrX'] 
    # chromosomes = ['chr22']

    special_value = np.Inf

    #for chromosome in chromosomes:
    hic_list = [process_chr(cell_type, args.chromosome, args.basedir, args.resolution, args.ref_scale, args.ref_gamma, special_value) for cell_type in cell_types]
    hic_list = [x for x in hic_list if x is not None]
    hic_list = [df.set_index(['bin1','bin2']) for df in hic_list]

    #Make average
    #Merge all hic matrices
    #Need to deal with nan vs 0 here. In the KR normalized matrices there are nan which we want to deal as missing. 
    #Rows that are not present in the hic dataframe should be considered 0
    #But after doing an outer join these rows will be represented as nan in the merged dataframe.
    #So need a way to distinguish nan vs 0.
    #Hack: convert all nan in the celltype specific hic dataframes to a special value. Then replace this special value after merging
    #TO DO: This is very memory intensive! (consider pandas.join or pandas.concat)
    # import pdb
    # pdb.set_trace()

    all_hic = pd.concat(hic_list, axis=1, join='outer', copy=False)
    hic_list = None #Clear from memory

    # import pdb
    # pdb.set_trace()

    #all_hic = pd.DataFrame().join(hic_list, how="outer", on=['bin1','bin2'])
    #all_hic = reduce(lambda x, y: pd.merge(x, y, on = ['bin1', 'bin2'], how = 'outer'), hic_list)
    all_hic.fillna(value=0, inplace=True)
    all_hic.replace(to_replace = special_value, value = np.nan, inplace=True)

    #compute the average
    cols_for_avg = list(filter(lambda x:'hic_kr' in x, all_hic.columns))
    
    #all_hic['avg_hic'] = all_hic[cols_for_avg].mean(axis=1)
    #avg_hic = all_hic[cols_for_avg].mean(axis=1)
    avg_hic = all_hic.mean(axis=1)
    num_good = len(cols_for_avg) - np.isnan(all_hic).sum(axis=1)

    #Check minimum number of cols
    all_hic.drop(cols_for_avg, inplace=True, axis=1)
    all_hic.reset_index(level=all_hic.index.names, inplace=True)
    all_hic['avg_hic'] = avg_hic.values
    all_hic.loc[num_good.values < args.min_cell_types_required, 'avg_hic'] = np.nan

    #Setup final matrix
    all_hic['bin1'] = all_hic['bin1'] * args.resolution
    all_hic['bin2'] = all_hic['bin2'] * args.resolution
    all_hic = all_hic.loc[np.logical_or(all_hic['avg_hic'] > 0, np.isnan(all_hic['avg_hic'])), ] # why do these 0's exist?

    os.makedirs(os.path.join(args.outDir, args.chromosome), exist_ok=True)
    all_hic.to_csv(os.path.join(args.outDir, args.chromosome, args.chromosome + ".avg.gz"), sep="\t", header=False, index=False, compression="gzip", na_rep=np.nan)

def scale_hic_with_powerlaw(hic, resolution, scale_ref, gamma_ref, scale, gamma):

    #get_powerlaw_at_distance expects positive gamma
    gamma_ref = -1 * gamma_ref
    gamma = -1 * gamma

    dists = (hic['bin2'] - hic['bin1']) * resolution
    pl_ref = get_powerlaw_at_distance(dists, gamma_ref, scale_ref) 
    pl = get_powerlaw_at_distance(dists, gamma, scale)

    hic['hic_kr'] = hic['hic_kr'] * (pl_ref / pl)

    return hic

def process_chr(cell_type, chromosome, basedir, resolution, scale_ref, gamma_ref, special_value):

    # import pdb
    # pdb.set_trace()

    hic_file, hic_norm_file, is_vc = get_hic_file(chromosome, os.path.join(basedir, cell_type, "5kb_resolution_intra"), allow_vc = True)
    if is_vc:
        return None
    # hic_file = os.path.join(basedir, cell_type, "5kb_resolution_intra", chromosome, chromosome + ".KRobserved")
    # hic_norm_file = os.path.join(basedir, cell_type, "5kb_resolution_intra", chromosome, chromosome + ".KRnorm")

    #Load gamma and scale
    pl_summary = pd.read_csv(os.path.join(basedir, cell_type, "5kb_resolution_intra/powerlaw/hic.powerlaw.txt"), sep="\t")

    #Read in and normalize to make DS
    hic = load_hic(hic_file = hic_file, 
                    hic_norm_file = hic_norm_file,
                    hic_is_vc=False,
                    hic_type="juicebox", 
                    hic_resolution = resolution, 
                    tss_hic_contribution=np.NaN, 
                    window = np.Inf, 
                    min_window=0, 
                    gamma = -1*pl_summary['pl_gamma'].values[0],
                    interpolate_nan=False, 
                    apply_diagonal_bin_correction=False)

    #power law scale 
    hic = scale_hic_with_powerlaw(hic, resolution, scale_ref, gamma_ref, scale = pl_summary['pl_scale'].values[0], gamma = pl_summary['pl_gamma'].values[0])

    #fill nan in hic matrix with special value
    #this will be turned back to nan after merging
    assert(not np.any(hic['hic_kr'].values == special_value))
    hic.loc[np.isnan(hic['hic_kr']), 'hic_kr'] = special_value

    return hic


if __name__ == '__main__':
    main()
