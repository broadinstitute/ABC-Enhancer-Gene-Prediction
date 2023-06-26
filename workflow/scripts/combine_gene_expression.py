import argparse
import pandas as pd
import numpy as np
import sys, traceback, os, os.path


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
        "--gene_expression_files",
        required=True,
        help="Genes to make predictions for. Formatted as the GeneList.txt file produced by run.neighborhoods.py",
    )
    parser.add_argument("--outdir", required=True, help="output directory")

    return parser


def get_predict_argument_parser():
    parser = get_model_argument_parser()
    return parser


def main():
    parser = get_predict_argument_parser()
    args = parser.parse_args()

    file_list = args.gene_expression_files.split(",")
    df_list = []
    for file in file_list:
        df = pd.read_csv(file, sep="\t", header=0)
        df_list.append(df)
    df_combined = pd.concat(df_list)
    print(df_combined.head(n=10))
    df_combined_wide = df_combined.pivot(
        columns="Sample",
        values="TargetGenePromoterActivityQuantile",
        index="TargetGene",
    ).add_prefix("TargetGenePromoterActivityQuantile_")
    df_combined_wide.to_csv(
        os.path.join(args.outdir, "Combined_gene_expression.txt"),
        sep="\t",
        index=True,
        header=True,
        float_format="%.6f",
        na_rep="NaN",
    )


if __name__ == "__main__":
    main()
