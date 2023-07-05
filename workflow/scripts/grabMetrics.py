#! bin/python3
import sys
from metrics import *
from tools import *
import glob
import argparse
import pickle
import os.path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--macs_peaks",
        help="narrowPeak file output by macs2. (eg. /users/kmualim/K562/macs2_peaks.narrowPeak)",
    )
    parser.add_argument("--preds_file", help="Prediction file output by ABC")
    parser.add_argument("--neighborhood_outdir", help="Neighborhood Directory")
    parser.add_argument("--outdir", help="Predictions Directory")
    args = parser.parse_args()
    return args


def generateQCMetrics(args):
    # prediction = "{}/EnhancerPredictionsFull.txt".format(args.preds_outdir)
    # read prediction file
    prediction_df = pd.read_csv(args.preds_file, sep="\t")
    # Generate QC Summary.txt in Predictions Directory
    pred_metrics = GrabQCMetrics(prediction_df, args.outdir)
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

    with open("{}/QCSummary.p".format(args.outdir), "wb") as f:
        pickle.dump(pred_metrics, f)


if __name__ == "__main__":
    args = parse_args()
    generateQCMetrics(args)
