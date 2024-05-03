import math
import time
from collections import defaultdict
from typing import Dict

import hicstraw  # hicstraw must be imported before pandas
import numpy as np
import pandas as pd
from hic import (
    get_hic_file,
    get_powerlaw_at_distance,
    load_hic_avg,
    load_hic_bedpe,
    load_hic_juicebox,
)
from tools import df_to_pyranges


def make_predictions(
    chromosome, enhancers, genes, args, hic_gamma, hic_scale, chrom_sizes_map
):
    pred = make_pred_table(chromosome, enhancers, genes, args.window, chrom_sizes_map)
    pred = annotate_predictions(pred, args.tss_slop)
    pred = add_powerlaw_to_predictions(pred, args, hic_gamma, hic_scale)
    # if Hi-C file is not provided, only powerlaw model will be computed
    if args.hic_file:
        if args.hic_type == "hic":
            pred = add_hic_from_hic_file(
                pred, args.hic_file, chromosome, args.hic_resolution
            )
        else:
            pred = add_hic_from_directory(
                chromosome,
                enhancers,
                genes,
                pred,
                args.hic_file,
                args,
                hic_gamma,
                hic_scale,
                chrom_sizes_map,
            )
        pred = qc_hic(pred, hic_gamma, hic_scale, args.hic_resolution)
        pred.fillna(value={"hic_contact": 0}, inplace=True)

        # Add powerlaw scaling
        pred = scale_hic_with_powerlaw(pred, args)
        # Add pseudocount
        pred = add_hic_pseudocount(pred, args)
        print("HiC Complete")

        pred = compute_score(
            pred,
            [pred["activity_base_enh"], pred["hic_contact_pl_scaled_adj"]],
            "ABC",
            adjust_self_promoters=True,
        )
    else:
        pred = compute_score(
            pred,
            [pred["activity_base_enh"], pred["powerlaw_contact"]],
            "ABC",
            adjust_self_promoters=True,
        )

    pred = compute_score(
        pred,
        [pred["activity_base_enh"], pred["powerlaw_contact"]],
        "powerlaw",
        adjust_self_promoters=True,
    )
    return pred


def make_pred_table(chromosome, enh, genes, window, chrom_sizes_map: Dict[str, int]):
    print("Making putative predictions table...")
    t = time.time()
    enh["enh_midpoint"] = (enh["start"] + enh["end"]) / 2
    enh["enh_idx"] = enh.index
    genes["gene_idx"] = genes.index
    enh_pr = df_to_pyranges(enh)
    genes_pr = df_to_pyranges(
        genes,
        start_col="TargetGeneTSS",
        end_col="TargetGeneTSS",
        start_slop=window,
        end_slop=window,
        chrom_sizes_map=chrom_sizes_map,
    )
    pred = enh_pr.join(genes_pr).df.drop(
        ["Start_b", "End_b", "chr_b", "Chromosome", "Start", "End"], axis=1
    )
    pred["distance"] = abs(pred["enh_midpoint"] - pred["TargetGeneTSS"])
    pred = pred.loc[pred["distance"] < window, :]  # for backwards compatability

    print(
        "Done. There are {} putative enhancers for chromosome {}".format(
            pred.shape[0], chromosome
        )
    )
    print("Elapsed time: {}".format(time.time() - t))

    return pred


def fill_diagonals(df, hic_resolution):
    """
    Fill diagonals based on the neighbors
    We want to search neighboring bins
    If hic_resolution is < 5kb, make sure we use a search space of 5kb for our
    neighbors
    """
    diagonal_bins = df[
        df.index.get_level_values("binX") == df.index.get_level_values("binY")
    ]
    search_space_bins = 1
    if hic_resolution < 5000:
        search_space_bins = math.ceil(5000 / hic_resolution)
    for (binX, binX), _ in diagonal_bins.iterrows():
        max_contact = 0
        for i in range(1, search_space_bins + 1):
            left_bin = (binX - i, binX)
            right_bin = (
                binX,
                binX + i,
            )  # we have to look above b/c we haven't processed right bin yet
            for bin in [left_bin, right_bin]:
                if bin in df.index:
                    max_contact = max(max_contact, df.loc[bin, "counts"])
        df.loc[(binX, binX), "counts"] = max_contact


def create_df_from_records(records, hic_resolution):
    """
    The bins returned from hic straw are in hic_resolution increments
    If resolution = 5k, we want to normalize the bins such that:
        5000 -> bin 1
        10000 -> bin 2
        ...
        100k -> bin 20
    Where the left number is the genomic position of the bin

    We also handle filling the diagonals differently
    """
    df = pd.DataFrame(records, columns=["binX", "binY", "counts"])
    df["binX"] = np.floor(df["binX"] / hic_resolution).astype(int)
    df["binY"] = np.floor(df["binY"] / hic_resolution).astype(int)
    df = df.set_index(["binX", "binY"])  # Set indexes for performance
    fill_diagonals(df, hic_resolution)
    return df


def get_chrom_format(hic: hicstraw.HiCFile, chromosome):
    """
    hic files can have 'chr1' or just '1' as the chromosome name
    we need to make sure we're using the format consistent with
    the hic file
    """
    hic_chrom_names = [chrom.name for chrom in hic.getChromosomes()]
    if hic_chrom_names[1].startswith("chr"):  # assume index 1 should be chr1
        return chromosome
    else:
        return chromosome[3:]


def get_chrom_size(hic: hicstraw.HiCFile, chromosome):
    for chrom in hic.getChromosomes():
        if chrom.name == chromosome:
            return chrom.length
    raise Exception(f"{chromosome} not found in hic data")


def determine_num_rows_to_fetch(
    chrom_size, hic_resolution, desired_mem_usage=16e9, record_size_bytes=200
):
    """
    Assuming each record takes up 200 bytes of memory, try to determine what is the
    max number of rows we can fetch from the hic matrix without going over the
    desired memory usage
    """
    num_records_per_row = math.ceil(chrom_size / hic_resolution)
    mem_usage_per_row = num_records_per_row * record_size_bytes
    max_rows = desired_mem_usage // mem_usage_per_row
    return int(max_rows)


def add_records_to_bin_sums(records, bin_sums, start, end):
    """
    hicstraw will remove duplicates in records
    e.g say there's contact between bin 1 and bin 2
    if i query for bin 1 and bin 2's rows in a single query, we will
    only get 1 value of (bin 1, bin 2).
    However, if we query for bin 1 and bin 2 separately, each result
    will give us (bin 1, bin 2).
    To correct for the duplicate values, we will divide a value in half
    if it'll get returned in a different query
    """
    for binX, binY, value in records:
        if binX < start or binY > end:
            value /= 2
        if binX == binY:
            bin_sums[binX] += value
        else:
            bin_sums[binX] += value
            bin_sums[binY] += value


def add_hic_from_hic_file(pred, hic_file, chromosome, hic_resolution):
    print("Begin HiC")
    start_time = time.time()
    pred["enh_bin"] = np.floor(pred["enh_midpoint"] / hic_resolution).astype(int)
    pred["tss_bin"] = np.floor(pred["TargetGeneTSS"] / hic_resolution).astype(int)
    pred["binX"] = np.min(pred[["enh_bin", "tss_bin"]], axis=1)
    pred["binY"] = np.max(pred[["enh_bin", "tss_bin"]], axis=1)
    hic = hicstraw.HiCFile(hic_file)
    chromosome = get_chrom_format(hic, chromosome)
    chromosome_size = get_chrom_size(hic, chromosome)
    matrix_object = hic.getMatrixZoomData(
        chromosome, chromosome, "observed", "SCALE", "BP", hic_resolution
    )
    # start and end loci need to cover the entire chromosome b/c we do normalization
    start_loci = 0
    end_loci = chromosome_size
    num_rows = 8000  # Using megamap with this number keeps memory usage right under 4GB
    bin_sums = defaultdict(float)
    step_size = num_rows * hic_resolution

    for i in range(start_loci, end_loci, step_size):
        start = i
        end = start + step_size - hic_resolution
        records = matrix_object.getRecords(start, end, start_loci, end_loci)
        if records:
            records = [[r.binX, r.binY, r.counts] for r in records]
            add_records_to_bin_sums(records, bin_sums, start, end)
            df = create_df_from_records(records, hic_resolution)
            pred = pred.merge(
                df,
                how="left",
                left_on=["binX", "binY"],
                right_index=True,
                suffixes=(None, "_"),
            )
            if "counts_" in pred:
                # After the first join, we have to merge the new records into
                # the existing count column. It's really just replacing the
                # nan values in pred_df with new values from latest records
                pred["counts"] = np.max(pred[["counts", "counts_"]], axis=1)
                pred.drop("counts_", inplace=True, axis=1)

    row_mean = np.mean(list(bin_sums.values()))
    # normalize hic_contact by the row_mean to reflect a doubly stochastic hic matrix
    pred["counts"] /= row_mean

    pred.drop(
        [
            "binX",
            "binY",
            "enh_idx",
            "gene_idx",
            "enh_midpoint",
            "tss_bin",
            "enh_bin",
            "max_neighbor_count",
        ],
        inplace=True,
        axis=1,
        errors="ignore",
    )
    print(
        "HiC added to predictions table. Elapsed time: {}".format(
            time.time() - start_time
        )
    )
    return pred.rename(columns={"counts": "hic_contact"})


def add_hic_from_directory(
    chromosome,
    enh,
    genes,
    pred,
    hic_dir,
    args,
    hic_gamma,
    hic_scale,
    chrom_sizes_map,
):
    hic_file, hic_norm_file, hic_is_vc = get_hic_file(
        chromosome, hic_dir, hic_type=args.hic_type
    )
    print("Begin HiC")
    # Add hic to pred table
    # At this point we have a table where each row is an enhancer/gene pair.
    # We need to add the corresponding HiC matrix entry.
    # If the HiC is provided in juicebox format (ie constant resolution), then we can just merge using the indices
    # But more generally we do not want to assume constant resolution. In this case hic should be provided in bedpe format

    t = time.time()
    if args.hic_type == "bedpe":
        HiC = load_hic_bedpe(hic_file)
        # Use pyranges to compute overlaps between enhancers/genes and hic bedpe table
        # Consider each range of the hic matrix separately - and merge each range into both enhancers and genes.
        # Then remerge on hic index

        HiC["hic_idx"] = HiC.index
        hic1 = df_to_pyranges(HiC, start_col="x1", end_col="x2", chr_col="chr1")
        hic2 = df_to_pyranges(HiC, start_col="y1", end_col="y2", chr_col="chr2")

        # Overlap in one direction
        enh_hic1 = (
            df_to_pyranges(
                enh,
                start_col="enh_midpoint",
                end_col="enh_midpoint",
                end_slop=1,
                chrom_sizes_map=chrom_sizes_map,
            )
            .join(hic1)
            .df
        )
        genes_hic2 = (
            df_to_pyranges(
                genes,
                start_col="TargetGeneTSS",
                end_col="TargetGeneTSS",
                end_slop=1,
                chrom_sizes_map=chrom_sizes_map,
            )
            .join(hic2)
            .df
        )
        ovl12 = enh_hic1[["enh_idx", "hic_idx", "hic_contact"]].merge(
            genes_hic2[["gene_idx", "hic_idx"]], on="hic_idx"
        )

        # Overlap in the other direction
        enh_hic2 = (
            df_to_pyranges(
                enh,
                start_col="enh_midpoint",
                end_col="enh_midpoint",
                end_slop=1,
                chrom_sizes_map=chrom_sizes_map,
            )
            .join(hic2)
            .df
        )
        genes_hic1 = (
            df_to_pyranges(
                genes,
                start_col="TargetGeneTSS",
                end_col="TargetGeneTSS",
                end_slop=1,
                chrom_sizes_map=chrom_sizes_map,
            )
            .join(hic1)
            .df
        )
        ovl21 = enh_hic2[["enh_idx", "hic_idx", "hic_contact"]].merge(
            genes_hic1[["gene_idx", "hic_idx"]], on=["hic_idx"]
        )

        # Concatenate both directions and merge into preditions
        ovl = pd.concat([ovl12, ovl21]).drop_duplicates()
        pred = pred.merge(ovl, on=["enh_idx", "gene_idx"], how="left")
    elif args.hic_type == "juicebox" or args.hic_type == "avg":
        if args.hic_type == "juicebox":
            HiC = load_hic_juicebox(
                hic_file=hic_file,
                hic_norm_file=hic_norm_file,
                hic_is_vc=hic_is_vc,
                hic_resolution=args.hic_resolution,
                tss_hic_contribution=args.tss_hic_contribution,
                window=args.window,
                min_window=0,
                gamma=hic_gamma,
                scale=hic_scale,
            )
        else:
            HiC = load_hic_avg(hic_file, args.hic_resolution)

        # Merge directly using indices
        # Could also do this by indexing into the sparse matrix (instead of merge) but this seems to be slower
        # Index into sparse matrix
        # pred['hic_contact'] = [HiC[i,j] for (i,j) in pred[['enh_bin','tss_bin']].values.tolist()]

        pred["enh_bin"] = np.floor(pred["enh_midpoint"] / args.hic_resolution).astype(
            int
        )
        pred["tss_bin"] = np.floor(pred["TargetGeneTSS"] / args.hic_resolution).astype(
            int
        )
        if not hic_is_vc:
            # in this case the matrix is upper triangular.
            #
            pred["bin1"] = np.amin(pred[["enh_bin", "tss_bin"]], axis=1)
            pred["bin2"] = np.amax(pred[["enh_bin", "tss_bin"]], axis=1)
            pred = pred.merge(HiC, how="left", on=["bin1", "bin2"])
        else:
            # The matrix is not triangular, its full
            # For VC assume genes correspond to rows and columns to enhancers
            pred = pred.merge(
                HiC,
                how="left",
                left_on=["tss_bin", "enh_bin"],
                right_on=["bin1", "bin2"],
            )

    pred.drop(
        [
            "x1",
            "x2",
            "y1",
            "y2",
            "bin1",
            "bin2",
            "enh_idx",
            "gene_idx",
            "hic_idx",
            "enh_midpoint",
            "tss_bin",
            "enh_bin",
        ],
        inplace=True,
        axis=1,
        errors="ignore",
    )

    print("HiC added to predictions table. Elapsed time: {}".format(time.time() - t))
    return pred


def scale_hic_with_powerlaw(pred, args):
    # Scale hic values to reference powerlaw

    if not args.scale_hic_using_powerlaw:
        #        values = pred.loc[pred['hic_contact']==0].index.astype('int')
        #        pred.loc[values, 'hic_contact'] = pred.loc[values, 'powerlaw_contact']
        pred["hic_contact_pl_scaled"] = pred["hic_contact"]
    else:
        pred["hic_contact_pl_scaled"] = pred["hic_contact"] * (
            pred["powerlaw_contact_reference"] / pred["powerlaw_contact"]
        )

    return pred


def add_powerlaw_to_predictions(pred, args, hic_gamma, hic_scale):
    pred["powerlaw_contact"] = get_powerlaw_at_distance(
        pred["distance"].values, hic_gamma, hic_scale
    )

    # 4.80 and 11.63 come from a linear regression of scale on gamma across 20
    # hic cell types at 5kb resolution. Do the params change across resolutions?
    hic_scale_reference = -4.80 + 11.63 * args.hic_gamma_reference
    pred["powerlaw_contact_reference"] = get_powerlaw_at_distance(
        pred["distance"].values, args.hic_gamma_reference, hic_scale_reference
    )
    return pred


def add_hic_pseudocount(pred, args):
    # Add a pseudocount based on the powerlaw expected count at a given distance
    pseudocount_distance = get_powerlaw_at_distance(
        args.hic_pseudocount_distance,
        args.hic_gamma,
        args.hic_scale,
        args.hic_pseudocount_distance,
    )
    pred["hic_pseudocount"] = pd.DataFrame(
        {"a": pred["powerlaw_contact"], "b": pseudocount_distance}
    ).min(axis=1)
    pred["hic_contact_pl_scaled_adj"] = (
        pred["hic_contact_pl_scaled"] + pred["hic_pseudocount"]
    )
    return pred


def qc_hic(pred, gamma, scale, resolution, threshold=0.01):
    # Gene promoters with insufficient hic coverage should get replaced with powerlaw
    summ = (
        pred.loc[pred["isSelfPromoter"], :]
        .groupby(["TargetGene"])
        .agg({"hic_contact": "sum"})
    )
    bad_genes = summ.loc[summ["hic_contact"] < threshold, :].index
    affected_pred = pred.loc[pred["TargetGene"].isin(bad_genes), :]

    pred.loc[affected_pred.index, "hic_contact"] = get_powerlaw_at_distance(
        affected_pred["distance"].values,
        gamma,
        scale,
        min_distance=resolution,
    )
    return pred


def compute_score(enhancers, product_terms, prefix, adjust_self_promoters=True):
    scores = np.column_stack(product_terms).prod(axis=1)

    enhancers[prefix + ".Score.Numerator"] = scores
    enhancers[prefix + ".Score"] = enhancers[
        prefix + ".Score.Numerator"
    ] / enhancers.groupby(["TargetGene", "TargetGeneTSS"])[
        prefix + ".Score.Numerator"
    ].transform(
        "sum"
    )

    # Self promoters by definition regulate the gene, so we
    # want to make sure they have a high score
    if adjust_self_promoters:
        self_promoters = enhancers[enhancers["isSelfPromoter"]]
        enhancers.loc[self_promoters.index, prefix + ".Score"] = 1

    return enhancers


def annotate_predictions(pred, tss_slop=500):
    # TO DO: Add is self genic
    pred["isSelfPromoter"] = np.logical_and.reduce(
        (
            pred["class"] == "promoter",
            pred.start - tss_slop < pred.TargetGeneTSS,
            pred.end + tss_slop > pred.TargetGeneTSS,
        )
    )

    return pred


def make_gene_prediction_stats(pred, score_column, threshold, output_file):
    summ1 = pred.groupby(["chr", "TargetGene", "TargetGeneTSS"]).agg(
        {
            "TargetGeneIsExpressed": lambda x: set(x).pop(),
            score_column: lambda x: all(np.isnan(x)),
            "name": "count",
        }
    )
    summ1.columns = ["geneIsExpressed", "geneFailed", "nEnhancersConsidered"]

    summ2 = (
        pred.loc[pred["class"] != "promoter", :]
        .groupby(["chr", "TargetGene", "TargetGeneTSS"])
        .agg({score_column: lambda x: sum(x > threshold)})
    )
    summ2.columns = ["nDistalEnhancersPredicted"]
    summ1 = summ1.merge(summ2, left_index=True, right_index=True)

    summ1.to_csv(output_file, sep="\t", index=True)
