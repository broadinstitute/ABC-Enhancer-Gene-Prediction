import argparse
import os
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all_putative")
    parser.add_argument("--score_column", default="ABC.Score")
    parser.add_argument("--chrom_sizes")
    parser.add_argument("--outdir")
    args = parser.parse_args()
    return args


def test_variant_overlap(args, all_putative):
    variant_overlap_file = os.path.join(
        args.outdir,
        "EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz",
    )
    # generate predictions for variant overlap
    score_t = all_putative[args.score_column] > 0.015
    not_promoter = all_putative["class"] != "promoter"
    is_promoter = all_putative["class"] == "promoter"
    score_one = all_putative[args.score_column] > 0.1
    all_putative[(score_t & not_promoter) | (is_promoter & score_one)]
    variant_overlap = all_putative[(score_t & not_promoter) | (is_promoter & score_one)]
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
        "zcat {}.tmp 2>/dev/null | head -1 | gzip > {}".format(
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
    all_putative = pd.read_csv(args.all_putative, sep="\t")
    test_variant_overlap(args, all_putative)
