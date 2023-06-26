import argparse
from predictor import *
from tools import *
from getVariantOverlap import *
import pandas as pd
import numpy as np
import sys, traceback, os, os.path
import time


def get_model_argument_parser():
    class formatter(
        argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
    ):
        pass

    parser = argparse.ArgumentParser(
        description="Predict enhancer relative effects.", formatter_class=formatter
    )
    readable = argparse.FileType("r")

    # Basic parameters
    parser.add_argument(
        "--genes",
        required=True,
        help="Genes to make predictions for. Formatted as the GeneList.txt file produced by run.neighborhoods.py",
    )
    parser.add_argument("--outdir", required=True, help="output directory")
    parser.add_argument(
        "--expression_cutoff",
        type=float,
        default=1,
        help="Make predictions for genes with expression higher than this value",
    )
    parser.add_argument(
        "--promoter_activity_quantile_cutoff",
        type=float,
        default=0,
        help="Quantile cutoff on promoter activity. Used to consider a gene 'expressed' in the absence of expression data",
    )
    parser.add_argument("--sampleid", required=True)

    return parser


def get_predict_argument_parser():
    parser = get_model_argument_parser()
    return parser


def main():
    parser = get_predict_argument_parser()
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    write_params(args, os.path.join(args.outdir, "parameters.predict.txt"))

    print("reading genes")
    genes = pd.read_csv(args.genes, sep="\t")
    genes = determine_expressed_genes(
        genes, args.expression_cutoff, args.promoter_activity_quantile_cutoff
    )
    genes = genes.loc[
        :,
        [
            "chr",
            "symbol",
            "tss",
            "Expression",
            "PromoterActivityQuantile",
            "isExpressed",
        ],
    ]
    genes.columns = [
        "chr",
        "TargetGene",
        "TargetGeneTSS",
        "TargetGeneExpression",
        "TargetGenePromoterActivityQuantile",
        "TargetGeneIsExpressed",
    ]

    genelist_slim = os.path.join(args.outdir, "GeneList_slim.txt")
    genelist_slim_expressed = os.path.join(args.outdir, "GeneListExpressed_slim.txt")
    genelist_slim_expressed_to_combine = os.path.join(
        args.outdir, "GeneListExpressed_slim_to_combine.txt"
    )
    print(genes.head(n=10))
    genes.to_csv(
        genelist_slim,
        sep="\t",
        index=False,
        header=True,
        float_format="%.6f",
        na_rep="NaN",
    )
    genes.loc[genes.TargetGeneIsExpressed, :].to_csv(
        genelist_slim_expressed,
        sep="\t",
        index=False,
        header=True,
        float_format="%.6f",
        na_rep="NaN",
    )
    genes_to_combine = genes.loc[
        genes.TargetGeneIsExpressed,
        ["TargetGene", "TargetGenePromoterActivityQuantile"],
    ]
    genes_to_combine["Sample"] = args.sampleid
    genes_to_combine.to_csv(
        genelist_slim_expressed_to_combine,
        sep="\t",
        index=False,
        header=True,
        float_format="%.6f",
        na_rep="NaN",
    )


if __name__ == "__main__":
    main()
