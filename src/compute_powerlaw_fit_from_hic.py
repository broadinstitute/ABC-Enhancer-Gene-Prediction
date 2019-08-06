#!/usr/bin/env python3

import numpy as np
import sys
import pandas
import scipy.sparse as ssp
import argparse
import glob
import os
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg')
import pylab

#To do: 
#1. Check if max/min window is off by 1 bin or is working properly
#2. Use MLE to compute powerlaw params and be done with scale parameter

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Helper to compute hic power-law fit parameters',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--bedDir', required=required_args, help="Directory containing bedgraphs. All files named *chr*.bg.gz will be loaded")
    parser.add_argument('--outDir', help="Output directory")
    parser.add_argument('--resolution', default=5000, help="Resolution of hic dataset (in bp)")
    parser.add_argument('--minWindow', default=10000, help="Minimum distance from gene TSS to compute normalizations (bp)")
    parser.add_argument('--maxWindow', default=1000000, help="Maximum distance from gene TSS to use to compute normalizations (bp)")

    args = parser.parse_args()
    return(args)

def welford(x):
    m = next(x)
    s = 0
    count = 1
    for k, val in enumerate(x):
        # k is 0 indexed, but we want k == natural index
        diff = val - m
        m_next = m + diff / (k + 2)
        s = s + diff.multiply(val - m_next)
        m = m_next
        count += 1

    return m, s / count


def filegen(args):
    maxsize = 5000000 + 10 * args.resolution 
    maxbin = 2 * maxsize // args.resolution + 1

    file_list = glob.glob(os.path.join(args.bedDir, "*chr*.bg.gz"))

    for f in file_list:
        try:
            df = pandas.read_table(f, compression='gzip', header=None)
        except pandas.io.common.EmptyDataError:
            # no data
            continue

        base = int(f.split('_')[-1].split('.')[0]) - maxsize

        # clip to bins
        locs = (df[1] - base) // args.resolution + 1
        valid = (locs < maxbin) & (locs >= 0)
        vals = df[3][valid]
        locs = locs[valid]
        norm = np.sum(vals)

        #data is all 0's
        if norm == 0:
            continue

        mat = ssp.csr_matrix((vals / norm, (0 * locs, locs)),
                             (1, 2 * maxsize // args.resolution + 1))
        yield mat


def compute_powerlaw_fit(m, args, make_plot=True):
    mean_hic = m.A[0]
    distance_from_center = abs(np.arange(len(mean_hic)) - len(mean_hic) // 2) + 1
    log_dist = np.log(distance_from_center)

    NBINS_AROUND_TSS_TO_REMOVE = int(args.resolution / args.minWindow)

    #Trim ends and remove regions close to the TSS
    midpoint = len(mean_hic)//2
    temp1 = np.arange(midpoint - args.maxWindow//args.resolution, midpoint - NBINS_AROUND_TSS_TO_REMOVE)
    temp2 = np.arange(midpoint + NBINS_AROUND_TSS_TO_REMOVE + 1, midpoint + args.maxWindow//args.resolution + 1)
    idx = np.concatenate((temp1,temp2))
    mean_hic_for_fit = mean_hic[idx]
    distance_from_center_for_fit = distance_from_center[idx]
    log_dist_for_fit = np.log(distance_from_center_for_fit)

    #Run linear regression in log log space
    def fit(params):
        gamma, scale_factor = params
        return np.log(mean_hic_for_fit) - (scale_factor + gamma * log_dist_for_fit)
    result = least_squares(fit, (-0.905, 0.25 * np.log(mean_hic_for_fit.max())))
    gamma, scale_factor = result.x

    # if make_plot:
    #     pylab.plot(np.log(mean_hic), 'k', label = 'data')
    #     pylab.plot(scale_factor + gamma * log_dist, 'r', label = 'fit')
    #     pylab.axvline(x=midpoint - args.maxWindow//args.resolution, color = 'b', linestyle='--', linewidth=1)
    #     pylab.axvline(x=midpoint + args.maxWindow//args.resolution, color = 'b', linestyle='--', linewidth=1)
    #     pylab.legend()
    #     #pylab.title(title_str)
    #     pylab.savefig(os.path.join(args.outDir, 'hic.powerlaw.png'))

    return result

if __name__ == '__main__':
    args = parseargs()

    os.makedirs(args.outDir, exist_ok=True)

    #Average together bedgraphs
    m, var = welford(filegen(args))

    #Save summary files
    np.savez(os.path.join(args.outDir, 'hic_bedgraph_summary.npz'), mean=m.A, var=var.A, resolution=args.resolution)
    pandas.DataFrame({ 'mean' : m.A[0], 'var' : var.A[0] }).to_csv(os.path.join(args.outDir, 'hic_bedgraph_summary.txt'), index=False, header=True, sep='\t')

    #compute normalization
    result = compute_powerlaw_fit(m, args)
    res = pandas.DataFrame({'resolution' : [args.resolution], 'maxWindow' : [args.maxWindow], 'minWindow' : [args.minWindow] ,'pl_gamma' : [result.x[0]], 'pl_scale' : [result.x[1]] })
    res.to_csv(os.path.join(args.outDir, 'hic.powerlaw.txt'), sep='\t', index=False, header=True)



