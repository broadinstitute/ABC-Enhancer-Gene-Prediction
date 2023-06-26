import argparse
import os
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--allPutative")
    parser.add_argument("--score_column", default="ABC.Score")
    parser.add_argument("--chrom_sizes")
    parser.add_argument("--outdir")
    parser.add_argument("--threshold")
    args = parser.parse_args()
    return args


def test_variant_overlap(args, all_putative):
    variant_overlap_file = os.path.join(
        args.outdir,
        f"EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp_ABC_threshold{args.threshold}.txt.gz",
    )
    # generate predictions for variant overlap
    score_t = all_putative[args.score_column] > float(args.threshold)
    not_promoter = all_putative["class"] != "promoter"
    variant_overlap = all_putative[(score_t & not_promoter)]
    # remove nan predictions
    variant_overlap_pred = variant_overlap.dropna(subset=[args.score_column])
    variant_overlap = variant_overlap_pred.loc[
        variant_overlap_pred["distance"] <= 2000000
    ]
    variant_overlap.to_csv(
        variant_overlap_file + ".tmp",
        sep="\t",
        index=False,
        header=True,
        compression="gzip",
        float_format="%.6f",
    )
    # shrink regions
    os.system(
        "zcat {}.tmp | head -1 | gzip > {}".format(
            variant_overlap_file, variant_overlap_file
        )
    )
    os.system(
        "zcat {}.tmp | sed 1d | bedtools slop -b -150 -g {} | gzip >> {}".format(
            variant_overlap_file, args.chrom_sizes, variant_overlap_file
        )
    )

    print("Done.")


if __name__ == "__main__":
    args = parse_args()
    all_putative = pd.read_csv(args.allPutative, sep="\t")
    test_variant_overlap(args, all_putative)
