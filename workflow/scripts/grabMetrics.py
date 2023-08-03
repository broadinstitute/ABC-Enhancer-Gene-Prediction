#! bin/python3
import argparse
import glob
import os.path
import pickle

import pandas as pd
from metrics import GrabQCMetrics, HiCQC, NeighborhoodFileQC, PeakFileQC


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--macs_peaks",
        required=True,
        help="narrowPeak file output by macs2. (eg. /users/kmualim/K562/macs2_peaks.narrowPeak)",
    )
    parser.add_argument(
        "--preds_file", required=True, help="Prediction file output by ABC"
    )
    parser.add_argument(
        "--neighborhood_outdir", required=True, help="Neighborhood Directory"
    )
    parser.add_argument("--chrom_sizes", required=True, help="Chromosome sizes file")
    parser.add_argument("--outdir", required=True, help="Predictions Directory")
    parser.add_argument(
        "--powerlaw_params_tsv",
        type=str,
        help="TSV file containing gamma/scale values according to powerlaw fit",
    )
    args = parser.parse_args()
    return args


def generateQCMetrics(args):
    chrom_order = pd.read_csv(args.chrom_sizes, sep="\t", header=None)[0].tolist()
    # read prediction file
    prediction_df = pd.read_csv(args.preds_file, sep="\t")
    # Generate QC Summary.txt in Predictions Directory
    pred_metrics = GrabQCMetrics(prediction_df, chrom_order, args.outdir)
    # Generate PeakFileQCSummary.txt in Peaks Directory#
    pred_metrics = PeakFileQC(pred_metrics, args.macs_peaks, args.outdir)
    # Appends Percentage Counts in Promoters into PeakFileQCSummary.txt
    if (
        len(
            glob.glob(
                os.path.join(
                    args.neighborhood_outdir, "Enhancers.DHS.*CountReads.bedgraph"
                )
            )
        )
        > 0
    ):
        pred_metrics = NeighborhoodFileQC(
            pred_metrics, args.neighborhood_outdir, args.outdir, "DHS"
        )
    if (
        len(
            glob.glob(
                os.path.join(
                    args.neighborhood_outdir, "Enhancers.H3K27ac.*CountReads.bedgraph"
                )
            )
        )
        > 0
    ):
        pred_metrics = NeighborhoodFileQC(
            pred_metrics, args.neighborhood_outdir, args.outdir, "H3K27ac"
        )
    if (
        len(
            glob.glob(
                os.path.join(
                    args.neighborhood_outdir, "Enhancers.ATAC.*CountReads.bedgraph"
                )
            )
        )
        > 0
    ):
        pred_metrics = NeighborhoodFileQC(
            pred_metrics, args.neighborhood_outdir, args.outdir, "ATAC"
        )
    if args.powerlaw_params_tsv:
        powerlaw_params = pd.read_csv(args.powerlaw_params_tsv, sep="\t").iloc[0]
        hic_gamma, hic_scale = (
            powerlaw_params["hic_gamma"],
            powerlaw_params["hic_scale"],
        )
        HiCQC(prediction_df, hic_gamma, hic_scale, args.outdir)

    with open("{}/QCSummary.p".format(args.outdir), "wb") as f:
        pickle.dump(pred_metrics, f)


if __name__ == "__main__":
    args = parse_args()
    generateQCMetrics(args)
