#! bin/python3
import argparse
import csv
import glob
import os.path

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
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
    parser.add_argument("--outdir", required=True, help="Metrics Directory")
    parser.add_argument("--output_qc_summary", required=True)
    parser.add_argument("--output_qc_plots", required=True)
    parser.add_argument(
        "--hic_gamma",
        type=float,
        help="Powerlaw exponent (gamma) to scale to. Must be positive",
    )
    parser.add_argument(
        "--hic_scale",
        type=float,
        help="scale of hic data. Must be positive",
    )
    args = parser.parse_args()
    return args


def generateQCMetrics(args):
    chrom_order = pd.read_csv(args.chrom_sizes, sep="\t", header=None)[0].tolist()
    # read prediction file
    prediction_df = pd.read_csv(args.preds_file, sep="\t")

    with PdfPages(args.output_qc_plots) as pdf_writer:
        pred_metrics = GrabQCMetrics(
            prediction_df, chrom_order, args.outdir, pdf_writer
        )
        pred_metrics = PeakFileQC(pred_metrics, args.macs_peaks, pdf_writer)
        if "hic_contact" in prediction_df.columns:
            HiCQC(prediction_df, args.hic_gamma, args.hic_scale, pdf_writer)

    # Appends Percentage Counts in Promoters into PeakFileQCSummary.txt
    potential_features = ["DHS", "H3K27ac", "ATAC"]
    for feature in potential_features:
        pattern = os.path.join(
            args.neighborhood_outdir, f"Enhancers.{feature}.*CountReads.bedgraph"
        )
        if glob.glob(pattern):
            pred_metrics = NeighborhoodFileQC(
                pred_metrics, args.neighborhood_outdir, feature
            )

    with open(args.output_qc_summary, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        for key, val in pred_metrics.items():
            writer.writerow((key, val))


if __name__ == "__main__":
    args = parse_args()
    generateQCMetrics(args)
