#!/usr/bin/env python3

import numpy as np
import pandas
import glob
import os.path
from intervaltree import IntervalTree
import pdb
import sys

class HiCFetcher(object):
    def __init__(self, dir, 
                 hic_gamma=1,
                 hic_gamma_reference=1,
                 scale_with_powerlaw=False, 
                 resolution=5000,
                 tss_hic_contribution=100):
        self.dir = dir
        self.hic_gamma = hic_gamma
        self.hic_gamma_reference = hic_gamma_reference
        self.resolution = resolution
        self.scale_with_powerlaw = scale_with_powerlaw
        self.needs_norm = True
        self.adjust_diag = True
        self.tss_hic_contribution = tss_hic_contribution

        # Fetch Hi-C data bedgraphs
        filenames = glob.glob(os.path.join(dir, '*chr*.bg.gz'))
        chroms = [f.split('.')[-3].split('_')[-2] for f in filenames]
        self._chromosomes = list(set(chroms))

        # Make intervaltree for each chromosome
        self.file_intervals = {ch: IntervalTree() for ch in chroms}
        for f, ch in zip(filenames, chroms):
            st = int(f.split('.')[-3].split('_')[-1])
            end = st + 1
            self.file_intervals[ch][st:end] = f

    def chromosomes(self):
        return self._chromosomes

    def query(self, chr, row, cols, enhancers, debug=False):
        intervals = list(self.file_intervals[chr][(row - self.resolution):(row + self.resolution)])
        if len(intervals) == 0:
            print("Could not find HiC data for {}:{}".format(chr, row))
            return np.full([len(cols), ], np.nan), np.nan, False, np.nan, np.nan

        if debug:
            print(intervals)

        # find best overlap
        best_interval = intervals[0]
        for i in intervals:
            if abs(row - (i.begin + i.end) / 2) < abs(row - (best_interval.begin + best_interval.end)):
                best_interval = i

        # load bedgraph    
        try:
            df = pandas.read_table(best_interval.data, compression='gzip', header=None)
            df.columns = ['chr', 'start', 'end', 'val']
        except:
            print("Could not load: " + best_interval.data)
            return np.full([len(cols), ], np.nan), np.nan, False, np.nan, np.nan
        
        #Adjust entry on the diagonal of the Hi-C matrix.
        if self.adjust_diag:
            diag_start = int(np.floor(best_interval[0] / self.resolution)*self.resolution)
            diag_end = int(np.ceil(best_interval[1] / self.resolution)*self.resolution)
            diag_idx = np.logical_and(df.start == diag_start, df.end == diag_end)
            assert(sum(diag_idx) <= 1)

            #Replace diagonal bin with max of neighboring bins multiplied by tss-scaling factor
            df.loc[diag_idx, 'val'] = df.loc[[diag_idx.idxmax() - 1, diag_idx.idxmax() + 1], 'val'].max() * self.tss_hic_contribution / 100

        #Normalize to sum to 1
        df.val /= df.val.sum()
        
        # find entries, handling missing data
        col_indices = np.searchsorted(df.start, cols, side='right') - 1
        valid = ((df.start[col_indices] <= cols) & (df.end[col_indices] > cols)).values
        values = np.zeros(len(cols))
        values[valid] = df.val[col_indices[valid]]
        rowmax = max(df.val)

        # Scale with respect to reference powerlaw 
        if self.scale_with_powerlaw:
            dists = abs(cols - row) / self.resolution
            log_dists = np.log(dists + 1)
            powerlaw_fit = -1*self.hic_gamma * log_dists
            powerlaw_fit_reference = -1*self.hic_gamma_reference * log_dists
            values_scaled = values * np.exp(powerlaw_fit_reference - powerlaw_fit)
            rowmax_scaled = rowmax
        else:
        	values_scaled = values
        	rowmax_scaled = rowmax

        return values_scaled, rowmax_scaled, True, values, rowmax

    def __call__(self, *args, **kwargs):
        return self.query(*args, **kwargs)

class DistanceModel(object):
    def __init__(self, model_gamma,
                 resolution=5000,
                 scale_with_powerlaw=True,
                 ):
        self.model_gamma = model_gamma
        self.resolution = resolution

    def __call__(self, *args, **kwargs):
        return self.query(*args, **kwargs)

    def query(self, dists):
        dists = (dists + 1) / self.resolution
        log_dists = np.log(dists + 1)
        return np.exp(-1*self.model_gamma * log_dists), 1

