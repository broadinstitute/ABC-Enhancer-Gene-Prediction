import argparse
import glob
import os
import sys
import traceback

import numpy as np
import pandas as pd
from hic import get_hic_file, load_hic_juicebox, load_hic_bedpe, load_hic_avg
from scipy import stats
from typing import Dict

# To do:
# 1. Use MLE to estimate exponent?
# 2. Support bedpe


def parseargs():
    class formatter(
        argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
    ):
        pass

    epilog = ""
    parser = argparse.ArgumentParser(
        description="Helper to compute hic power-law fit parameters",
        epilog=epilog,
        formatter_class=formatter,
    )
    readable = argparse.FileType("r")
    parser.add_argument(
        "--hic_dir",
        help="Directory containing observed HiC KR normalized matrices. File naming and structure should be: hicDir/chr*/chr*.{KR,Interscale}observed",
    )
    parser.add_argument("--outDir", help="Output directory")
    parser.add_argument(
        "--hic_type",
        default="juicebox",
        choices=["juicebox", "bedpe", "avg"],
        help="format of hic files",
    )
    parser.add_argument(
        "--hic_resolution",
        default=5000,
        type=int,
        help="For Juicebox: resolution of hic dataset (in bp). For bedpe: distances will be binned to this resolution for powerlaw fit",
    )
    parser.add_argument(
        "--minWindow",
        default=5000,
        type=int,
        help="Minimum distance between bins to include in powerlaw fit (bp). Recommended to be at least >= resolution to avoid using the diagonal of the HiC Matrix",
    )
    parser.add_argument(
        "--maxWindow",
        default=1000000,  # 1Mbp
        type=int,
        help="Maximum distance between bins to include in powerlaw fit (bp)",
    )
    parser.add_argument(
        "--chr",
        default="all",
        help="Comma delimited list of chromosomes to use for fit. Defualts to chr[1..22],chrX",
    )

    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    os.makedirs(args.outDir, exist_ok=True)

    if args.chr == "all":
        chromosomes = ["chr" + str(x) for x in list(range(1, 23))] + ["chrX"]
    else:
        chromosomes = args.chr.split(",")

    HiC = load_hic_for_powerlaw(
        chromosomes,
        args.hic_dir,
        args.hic_type,
        args.hic_resolution,
        args.minWindow,
        args.maxWindow,
    )

    # Run
    slope, intercept, hic_mean_var = do_powerlaw_fit(HiC, args.hic_resolution)

    # print
    res = pd.DataFrame(
        {
            "resolution": [args.hic_resolution],
            "maxWindow": [args.maxWindow],
            "minWindow": [args.minWindow],
            "hic_gamma": [slope * -1],  # gamma defined as neg slope
            "hic_scale": [intercept],
        }
    )
    res.to_csv(
        os.path.join(args.outDir, "hic.powerlaw.tsv"),
        sep="\t",
        index=False,
        header=True,
    )

    hic_mean_var.to_csv(
        os.path.join(args.outDir, "hic.mean_var.tsv"), sep="\t", index=True, header=True
    )


def load_hic_for_powerlaw(
    chromosomes, hicDir, hic_type, hic_resolution, min_window, max_window
):
    all_data_list = []
    for chrom in chromosomes:
        try:
            hic_file, hic_norm_file, hic_is_vc = get_hic_file(
                chrom, hicDir, hic_type=hic_type, allow_vc=False
            )
            if hic_type == "juicebox":
                print("Working on {}".format(hic_file))
                this_data = load_hic_juicebox(
                    hic_file=hic_file,
                    hic_norm_file=hic_norm_file,
                    hic_is_vc=hic_is_vc,
                    hic_resolution=hic_resolution,
                    tss_hic_contribution=100,
                    window=max_window,
                    min_window=min_window,
                    gamma=np.nan,
                    scale=np.nan,
                    interpolate_nan=False,
                )
                this_data["dist_for_fit"] = (
                    abs(this_data["bin1"] - this_data["bin2"]) * hic_resolution
                )
            elif hic_type == "bedpe":
                print("Working on {}".format(hic_file))
                this_data = load_hic_bedpe(hic_file)
                # Compute distance in bins as with juicebox data.
                # This is needed to in order to maintain consistency, but is probably slightly less accurate.
                # Binning also reduces noise level.
                rawdist = abs(
                    (this_data["x2"] + this_data["x1"]) / 2
                    - (this_data["y2"] + this_data["y1"]) / 2
                )
                this_data["dist_for_fit"] = (rawdist // hic_resolution) * hic_resolution
                this_data = this_data.loc[
                    np.logical_and(
                        this_data["dist_for_fit"] >= min_window,
                        this_data["dist_for_fit"] <= max_window,
                    )
                ]
            elif hic_type == "avg":
                print("Working on {}".format(hic_file))
                this_data = load_hic_avg(hic_file, hic_resolution)
                this_data["dist_for_fit"] = (
                    abs((this_data["bin1"] - this_data["bin2"]) / hic_resolution)
                    * hic_resolution
                )
                this_data["dist_for_fit"] = this_data["dist_for_fit"].astype("int")
            else:
                raise Exception("invalid --hic_type")

            all_data_list.append(this_data)
        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)

    all_data = pd.concat(all_data_list)

    return all_data


def do_powerlaw_fit(HiC, resolution):
    print("Running regression")

    # TO DO:
    # Print out mean/var plot of powerlaw relationship
    # Juicebox output is in sparse matrix format. This is an attempt to get a "mean" hic_contact for each distance bin since we don't have the total # of bins available.
    HiC_summary = HiC.groupby("dist_for_fit").agg({"hic_contact": "sum"})
    # HiC_summary['hic_contact'] = HiC_summary.hic_contact / HiC_summary.hic_contact.sum() #technically this normalization should be over the entire genome (not just to maxWindow). Will only affect intercept though
    # get a better approximation of the total # of bins .
    total_num_bins = len(HiC.loc[HiC["dist_for_fit"] == resolution])
    print(total_num_bins)
    #    idx = np.isfinite(HiC_summary['hic_contact']) & np.isfinite(HiC_summary.index)
    HiC_summary["hic_contact"] = HiC_summary.hic_contact / total_num_bins
    pseudocount = 0.000001
    #    pseudocount = 0
    res = stats.linregress(
        np.log(HiC_summary.index + pseudocount),
        np.log(HiC_summary["hic_contact"] + pseudocount),
    )

    hic_mean_var = HiC.groupby("dist_for_fit").agg({"hic_contact": ["mean", "var"]})
    hic_mean_var.columns = ["mean", "var"]

    return res.slope, res.intercept, hic_mean_var


if __name__ == "__main__":
    main()
