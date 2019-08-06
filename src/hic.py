#!/usr/bin/env python3

from weakref import WeakValueDictionary

import numpy as np
import scipy.sparse as ssp
import pandas

class TempDict(dict):
    pass

class HiC(object):
    def __init__(self, files, window=5000000, resolution=5000, kr_cutoff=0.1):
        self.files = files
        self.cache = WeakValueDictionary()
        self.window = window
        self.resolution = resolution
        self.kr_cutoff = kr_cutoff
        self.__last = None  # to keep reference to last matrix (avoid kicking out of cache)

        self._chromosomes = list(files.keys())
        assert len(self._chromosomes) > 0, "No HiC data found"

    def chromosomes(self):
        return self._chromosomes

    def row(self, chr, row):
        try:
            hic = self.cache[chr]
        except KeyError:
            hic = self.load(chr)
            self.__last = self.cache[chr] = hic

        hicdata = hic['hic_mat']
        norms = hic['hic_norm']

        # find row in matrix
        rowidx = row // self.resolution

        # clip to range for which we have data
        rowidx = max(0, min(rowidx, hicdata.shape[0] - 1))

        #Set all entries that have nan normalization factor to nan 
        data = hicdata[rowidx, :]
        if norms is not None:
            data[0, np.where(np.isnan(norms))] = np.nan

        return data

    def query(self, chr, row, cols):
        try:
            hicdata = self.cache[chr]
        except KeyError:
            hicdata = self.load(chr)
            self.__last = self.cache[chr] = hicdata

        # find cols in matrix
        colsidx = cols // self.resolution
        valid_colsidx = np.clip(colsidx, 0, hicdata.shape[1] - 1)
        rowdata = self.row(row)

        # extract column values
        values = rowdata[:, valid_colsidx].todense().A.ravel()

        # out-of-bound values == 0
        values[colsidx != valid_colsidx] = 0

        return values

    def __call__(self, *args, **kwargs):
        return self.query(*args, **kwargs)

    def load(self, chr):

        #TO DO: What is this?
        hic_filename = self.files[chr]
        norm_filename = None
        if isinstance(hic_filename, tuple):
            hic_filename, norm_filename = hic_filename

        print("loading", hic_filename)
        sparse_matrix = hic_to_sparse(hic_filename,
                                      self.window, self.resolution)

        if norm_filename is not None:
            norms = np.loadtxt(norm_filename)
            assert len(norms) >= sparse_matrix.shape[0]
            if len(norms) > sparse_matrix.shape[0]:
                norms = norms[:sparse_matrix.shape[0]] #JN: 4/23/18 - is this always guaranteed to be correct???

            norms[norms < self.kr_cutoff] = np.nan
            norm_mat = ssp.dia_matrix((1.0 / norms, [0]), (len(norms), len(norms)))

            # normalize row and columns
            sparse_matrix_norm = norm_mat * sparse_matrix * norm_mat

        return TempDict(hic_mat=sparse_matrix_norm, hic_norm=norms)


def hic_to_sparse(filename, window, resolution):
    HiC = pandas.read_table(filename, names=["start", "end", "counts"],
                            header=None, engine='c', memory_map=True)

    # verify our assumptions
    assert np.all(HiC.start <= HiC.end)

    # find largest entry
    max_pos = max(HiC.start.max(), HiC.end.max())
    hic_size = max_pos // resolution + 1

    # drop NaNs from hic
    HiC = HiC[~np.isnan(HiC.counts.as_matrix())]
    print("HiC has {} rows after dropping NaNs".format(HiC.shape[0]))

    # window distance between contacts
    too_far = (HiC.end - HiC.start) >= window
    HiC = HiC[~too_far]
    print("HiC has {} rows after windowing to {}".format(HiC.shape[0], window))

    # convert to sparse matrix in CSR (compressed sparse row) format, chopping
    # down to HiC bin size.  note that conversion to scipy sparse matrices
    # accumulates repeated indices, so this will do the right thing.
    row = np.floor(HiC.start.as_matrix() / resolution).astype(int)
    col = np.floor(HiC.end.as_matrix() / resolution).astype(int)
    dat = HiC.counts.as_matrix()
    # we want a symmetric matrix.  Easiest to do that during creation, but have to be careful of diagonal
    mask = (row != col)  # off-diagonal
    row2 = col[mask]  # note the row/col swap
    col2 = row[mask]
    dat2 = dat[mask]

    # concat and create
    row = np.hstack((row, row2))
    col = np.hstack((col, col2))
    dat = np.hstack((dat, dat2))
    return ssp.csr_matrix((dat, (row, col)), (hic_size, hic_size))
