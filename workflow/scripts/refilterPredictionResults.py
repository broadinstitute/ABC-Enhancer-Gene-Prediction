import argparse
from predictor import *
from tools import *
from getVariantOverlap import *
import pandas as pd
import numpy as np
import sys, traceback, os, os.path
import time


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--allPutative", help="Peaks Directory")
    parser.add_argument("--alternative_threshold", help="Neighborhood Directory")
    parser.add_argument("--score_column", default="ABC.Score"),
    parser.add_argument("--outdir", help="Output Directory")
    args = parser.parse_args()
    return args


def filter(args):
    pred_file_full = os.path.join(
        args.outdir,
        f"Predictions_ABC_threshold{args.alternative_threshold}",
        "EnhancerPredictionsFull.txt",
    )
    pred_file_bedpe = os.path.join(
        args.outdir,
        f"Predictions_ABC_threshold{args.alternative_threshold}",
        "EnhancerPredictions.bedpe",
    )
    all_putative_path = args.allPutative
    all_putative = pd.read_csv(all_putative_path, compression="gzip", sep="\t")
    all_positive = all_putative.iloc[
        np.logical_and.reduce(
            (
                all_putative[args.score_column] > float(args.alternative_threshold),
                ~(all_putative["class"] == "promoter"),
            )
        ),
        :,
    ]
    all_positive.to_csv(
        pred_file_full, sep="\t", index=False, header=True, float_format="%.6f"
    )
    write_connections_bedpe_format(all_positive, pred_file_bedpe, args.score_column)


if __name__ == "__main__":
    args = parse_args()
    filter(args)
