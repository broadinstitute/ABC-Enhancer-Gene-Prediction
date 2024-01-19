import os
import time

import numpy as np
import pandas as pd
import scipy.sparse as ssp


def get_hic_file(chromosome, hic_dir, allow_vc=True, hic_type="juicebox"):
    if hic_type == "juicebox":
        is_vc = False
        filetypes = ["KR", "INTERSCALE"]
        for filetype in filetypes:
            hic_file = os.path.join(
                hic_dir, chromosome, chromosome + f".{filetype}observed.gz"
            )
            hic_norm = os.path.join(
                hic_dir, chromosome, chromosome + f".{filetype}norm.gz"
            )
            if hic_exists(hic_file):
                print("Using: " + hic_file)
                return hic_file, hic_norm, False

        if allow_vc:
            hic_file = os.path.join(hic_dir, chromosome, chromosome + ".VCobserved.gz")
            hic_norm = os.path.join(hic_dir, chromosome, chromosome + ".VCnorm.gz")
            if hic_exists(hic_file):
                print(
                    f"Could not find KR normalized hic file. Using VC normalized hic file: {hic_file}"
                )
                return hic_file, hic_norm, True

        raise RuntimeError(
            f"Could not find {', '.join(filetypes)} or VC normalized hic files"
        )

    elif hic_type == "bedpe":
        hic_file = os.path.join(hic_dir, chromosome, chromosome + ".bedpe.gz")
        return hic_file, None, None
    elif hic_type == "avg":
        hic_file = os.path.join(hic_dir, chromosome, chromosome + ".bed.gz")
        return hic_file, None, None


def hic_exists(file):
    if not os.path.exists(file):
        return False
    elif file.endswith("gz"):
        # gzip file still have some size. This is a hack
        return os.path.getsize(file) > 100
    else:
        return os.path.getsize(file) > 0


def load_hic_juicebox(
    hic_file,
    hic_norm_file,
    hic_is_vc,
    hic_resolution,
    tss_hic_contribution,
    window,
    min_window,
    gamma,
    scale=None,
    apply_diagonal_bin_correction=True,
):
    print("Loading HiC Juicebox")
    HiC_sparse_mat = hic_to_sparse(hic_file, hic_norm_file, hic_resolution)
    HiC = process_hic(
        hic_mat=HiC_sparse_mat,
        hic_norm_file=hic_norm_file,
        hic_is_vc=hic_is_vc,
        resolution=hic_resolution,
        tss_hic_contribution=tss_hic_contribution,
        window=window,
        min_window=min_window,
        gamma=gamma,
        apply_diagonal_bin_correction=apply_diagonal_bin_correction,
        scale=scale,
    )
    return HiC


def load_hic_bedpe(hic_file):
    print("Loading HiC bedpe")
    return pd.read_csv(
        hic_file,
        sep="\t",
        names=["chr1", "x1", "x2", "chr2", "y1", "y2", "name", "hic_contact"],
    )


def load_hic_avg(hic_file, hic_resolution):
    print("Loading HiC avg")
    cols = {"x1": np.int64, "x2": np.int64, "hic_contact": np.float64}
    HiC = pd.read_csv(
        hic_file, sep="\t", names=cols.keys(), usecols=cols.keys(), dtype=cols
    )
    HiC["x1"] = np.floor(HiC["x1"] / hic_resolution).astype(int)
    HiC["x2"] = np.floor(HiC["x2"] / hic_resolution).astype(int)
    HiC.rename(columns={"x1": "bin1", "x2": "bin2"}, inplace=True)
    return HiC


# def juicebox_to_bedpe(hic, chromosome, resolution):
#     hic['chr'] = chromosome
#     hic['x1'] = hic['bin1'] * resolution
#     hic['x2'] = (hic['bin1'] + 1) * resolution
#     hic['y1'] = hic['bin2'] * resolution
#     hic['y2'] = (hic['bin2'] + 1) * resolution

#     return(hic)


def process_hic(
    hic_mat,
    hic_norm_file,
    hic_is_vc,
    resolution,
    tss_hic_contribution,
    window,
    min_window=0,
    hic_is_doubly_stochastic=False,
    apply_diagonal_bin_correction=True,
    gamma=None,
    kr_cutoff=0.25,
    scale=None,
):
    # Make doubly stochastic.
    # Juicer produces a matrix with constant row/column sums. But sum is not 1 and is variable across chromosomes
    t = time.time()

    if not hic_is_doubly_stochastic and not hic_is_vc:
        # Any row with Nan in it will sum to nan
        # So need to calculate sum excluding nan
        temp = hic_mat
        temp.data = np.nan_to_num(temp.data, copy=False)
        sums = temp.sum(axis=0)
        sums = sums[~np.isnan(sums)]
        # assert(np.max(sums[sums > 0])/np.min(sums[sums > 0]) < 1.001)
        mean_sum = np.mean(sums[sums > 0])

        if abs(mean_sum - 1) < 0.001:
            print(
                "HiC Matrix has row sums of {}, continuing without making doubly stochastic".format(
                    mean_sum
                )
            )
        else:
            print(
                "HiC Matrix has row sums of {}, making doubly stochastic...".format(
                    mean_sum
                )
            )
            hic_mat = hic_mat.multiply(1 / mean_sum)

    # Adjust diagonal of matrix based on neighboring bins
    # First and last rows need to be treated differently
    if apply_diagonal_bin_correction:
        last_idx = hic_mat.shape[0] - 1
        nonzero_diag = hic_mat.nonzero()[0][
            hic_mat.nonzero()[0] == hic_mat.nonzero()[1]
        ]
        nonzero_diag = list(
            set(nonzero_diag) - set(np.array([last_idx])) - set(np.array([0]))
        )

        for ii in nonzero_diag:
            hic_mat[ii, ii] = (
                max(hic_mat[ii, ii - 1], hic_mat[ii, ii + 1])
                * tss_hic_contribution
                / 100
            )

        if hic_mat[0, 0] != 0:
            hic_mat[0, 0] = hic_mat[0, 1] * tss_hic_contribution / 100

        if hic_mat[last_idx, last_idx] != 0:
            hic_mat[last_idx, last_idx] = (
                hic_mat[last_idx, last_idx - 1] * tss_hic_contribution / 100
            )

    # Remove lower triangle
    if not hic_is_vc:
        hic_mat = ssp.triu(hic_mat)
    else:
        hic_mat = process_vc(hic_mat)

    # Turn into dataframe
    hic_mat = hic_mat.tocoo(copy=False)
    hic_df = pd.DataFrame(
        {"bin1": hic_mat.row, "bin2": hic_mat.col, "hic_contact": hic_mat.data}
    )

    # Prune to window
    hic_df = hic_df.loc[
        np.logical_and(
            abs(hic_df["bin1"] - hic_df["bin2"]) <= window / resolution,
            abs(hic_df["bin1"] - hic_df["bin2"]) >= min_window / resolution,
        )
    ]
    print(
        "HiC has {} rows after windowing between {} and {}".format(
            hic_df.shape[0], min_window, window
        )
    )

    print("process.hic: Elapsed time: {}".format(time.time() - t))

    return hic_df


def hic_to_sparse(filename, norm_file, resolution, hic_is_doubly_stochastic=False):
    t = time.time()
    HiC = pd.read_table(
        filename,
        names=["bin1", "bin2", "hic_contact"],
        header=None,
        engine="c",
        memory_map=True,
    )

    # verify our assumptions
    assert np.all(HiC.bin1 <= HiC.bin2)

    # Need load norms here to know the dimensions of the hic matrix
    norms = pd.read_csv(norm_file, header=None)
    hic_size = norms.shape[0]

    # convert to sparse matrix in CSR (compressed sparse row) format, chopping
    # down to HiC bin size.  note that conversion to scipy sparse matrices
    # accumulates repeated indices, so this will do the right thing.
    row = np.floor(HiC.bin1.values / resolution).astype(int)
    col = np.floor(HiC.bin2.values / resolution).astype(int)
    dat = HiC.hic_contact.values

    # JN: Need both triangles in order to compute row/column sums to make double stochastic.
    # If juicebox is upgraded to return DS matrices, then can remove one triangle
    # TO DO: Remove one triangle when juicebox is updated.
    # we want a symmetric matrix.  Easiest to do that during creation, but have to be careful of diagonal
    if not hic_is_doubly_stochastic:
        mask = row != col  # off-diagonal
        row2 = col[mask]  # note the row/col swap
        col2 = row[mask]
        dat2 = dat[mask]

        # concat and create
        row = np.hstack((row, row2))
        col = np.hstack((col, col2))
        dat = np.hstack((dat, dat2))

    print("hic.to.sparse: Elapsed time: {}".format(time.time() - t))

    return ssp.csr_matrix((dat, (row, col)), (hic_size, hic_size))


def get_powerlaw_at_distance(distances, gamma, scale, min_distance=5000):
    assert gamma > 0
    assert scale > 0

    # The powerlaw is computed for distances > 5kb. We don't know what the contact freq looks like at < 5kb.
    # So just assume that everything at < 5kb is equal to 5kb.
    # TO DO: get more accurate powerlaw at < 5kb
    distances = np.clip(distances, min_distance, np.Inf)
    log_dists = np.log(distances + 1)

    powerlaw_contact = np.exp(scale + -1 * gamma * log_dists)
    return powerlaw_contact


def process_vc(hic):
    # For a vc normalized matrix, need to make rows sum to 1.
    # Assume rows correspond to genes and cols to enhancers

    row_sums = hic.sum(axis=0)
    row_sums[row_sums == 0] = 1
    norm_mat = ssp.dia_matrix(
        (1.0 / row_sums, [0]), (row_sums.shape[1], row_sums.shape[1])
    )

    # left multiply to operate on rows
    hic = norm_mat * hic

    return hic
