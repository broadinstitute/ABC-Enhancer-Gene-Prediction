import argparse
import os
import os.path
import sys
import time
import traceback
from typing import Dict

import numpy as np
import pandas as pd
from compute_powerlaw_fit_from_hic import do_powerlaw_fit, load_hic_for_powerlaw
from getVariantOverlap import *
from predictor import *
from tools import *


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
        "--enhancers",
        required=True,
        help="Candidate enhancer regions. Formatted as the EnhancerList.txt file produced by run.neighborhoods.py",
    )
    parser.add_argument(
        "--genes",
        required=True,
        help="Genes to make predictions for. Formatted as the GeneList.txt file produced by run.neighborhoods.py",
    )
    parser.add_argument("--outdir", required=True, help="output directory")
    parser.add_argument(
        "--score_column",
        default="ABC.Score",
        help="Column name of score to use for thresholding",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        required=True,
        default=0.022,
        help="Threshold on ABC Score (--score_column) to call a predicted positive. Note that the threshold will need to be adjusted based on the combination of input datasets used.",
    )
    parser.add_argument("--cellType", help="Name of cell type")
    parser.add_argument("--chrom_sizes", required=True, help="Chromosome sizes file")

    # hic
    # To do: validate params
    parser.add_argument("--hic_dir", default=None, help="HiC directory")
    parser.add_argument("--hic_resolution", type=int, help="HiC resolution")
    parser.add_argument(
        "--hic_pseudocount_distance",
        type=int,
        default=1e6,
        help="A pseudocount is added equal to the powerlaw fit at this distance",
    )
    parser.add_argument(
        "--hic_type",
        default="juicebox",
        choices=["juicebox", "bedpe", "avg"],
        help="format of hic files",
    )
    parser.add_argument(
        "--hic_is_doubly_stochastic",
        action="store_true",
        help="If hic matrix is already doubly stochastic, can skip this step",
    )

    # Power law
    parser.add_argument(
        "--scale_hic_using_powerlaw",
        action="store_true",
        help="Quantile normalize Hi-C values using powerlaw relationship. This parameter will rescale Hi-C contacts from the input Hi-C data (specified by --hic_gamma and --hic_scale) to match the power-law relationship of a reference cell type (specified by --hic_gamma_reference)",
    )
    parser.add_argument(
        "--powerlaw_params_tsv",
        type=str,
        help="TSV file containing gamma/scale values according to powerlaw fit",
    )
    parser.add_argument(
        "--hic_gamma_reference",
        type=float,
        default=0.87,
        help="Powerlaw exponent (gamma) to scale to. Must be positive",
    )

    # Genes to run through model
    parser.add_argument(
        "--run_all_genes",
        action="store_true",
        help="Do not check for gene expression, make predictions for all genes. Use of this parameter is not recommended.",
    )
    parser.add_argument(
        "--expression_cutoff",
        type=float,
        default=1,
        help="Make predictions for genes with expression higher than this value. Use of this parameter is not recommended.",
    )
    parser.add_argument(
        "--promoter_activity_quantile_cutoff",
        type=float,
        default=0.4,
        help="Quantile cutoff on promoter activity. Used to consider a gene 'expressed' in the absence of expression data",
    )

    # Output formatting
    parser.add_argument(
        "--make_all_putative",
        action="store_true",
        help="Make big file with concatenation of all genes file",
    )
    parser.add_argument(
        "--use_hdf5",
        action="store_true",
        help="Write AllPutative file in hdf5 format instead of tab-delimited",
    )

    # Parameters used in development of ABC model that should not be changed for most use cases
    parser.add_argument(
        "--window",
        type=int,
        default=5000000,
        help="Consider all candidate elements within this distance of the gene's TSS when computing ABC scores. This was a parameter optimized during development of ABC and should not typically be changed.",
    )
    parser.add_argument(
        "--tss_hic_contribution",
        type=float,
        default=100,
        help="Diagonal bin of Hi-C matrix is set to this percentage of the maximum of its neighboring bins. Default value (100%) means that the diagonal bin of the Hi-C matrix will be set exactly to the maximum of its two neighboring bins. This is a parameter used in development of the ABC model that should not typically be changed, unless experimenting with using different types of input 3D contact data that require different handling of the diagonal bin.",
    )

    # Other
    parser.add_argument(
        "--tss_slop",
        type=int,
        default=500,
        help="Distance from tss to search for self-promoters",
    )
    parser.add_argument(
        "--chromosomes",
        default="all",
        help="chromosomes to make predictions for. Defaults to intersection of all chromosomes in --genes and --enhancers",
    )
    parser.add_argument(
        "--include_chrY",
        "-y",
        action="store_true",
        help="Make predictions on Y chromosome",
    )
    parser.add_argument(
        "--include_self_promoter",
        action="store_true",
        help="Include self-promoter elements as enhancers ",
    )

    return parser


def get_predict_argument_parser():
    parser = get_model_argument_parser()
    return parser


def main():
    parser = get_predict_argument_parser()
    args = parser.parse_args()

    validate_args(args)

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
            "DHS.RPKM.quantile.TSS1Kb",
        ],
    ]
    genes.columns = [
        "chr",
        "TargetGene",
        "TargetGeneTSS",
        "TargetGeneExpression",
        "TargetGenePromoterActivityQuantile",
        "TargetGeneIsExpressed",
        "normalized_dhs",
    ]

    print("reading enhancers")
    enhancers_full = pd.read_csv(args.enhancers, sep="\t")
    # TO DO
    # Think about which columns to include
    enhancers = enhancers_full.loc[
        :,
        [
            "chr",
            "start",
            "end",
            "name",
            "class",
            "activity_base",
            "normalized_dhs",
            "normalized_h3K27ac",
        ],
    ]
    enhancers["activity_base_squared"] = enhancers["activity_base"] ** 2
    # Initialize Prediction files
    pred_file_full = os.path.join(args.outdir, "EnhancerPredictionsFull.tsv")
    pred_file_slim = os.path.join(args.outdir, "EnhancerPredictions.tsv")
    pred_file_bedpe = os.path.join(args.outdir, "EnhancerPredictions.bedpe")
    all_pred_file_expressed = os.path.join(
        args.outdir, "EnhancerPredictionsAllPutative.txt.gz"
    )
    all_pred_file_nonexpressed = os.path.join(
        args.outdir, "EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz"
    )
    variant_overlap_file = os.path.join(
        args.outdir,
        "EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.txt.gz",
    )
    all_putative_list = []

    # Make predictions
    if args.chromosomes == "all":
        chromosomes = set(genes["chr"]).intersection(set(enhancers["chr"]))
        if not args.include_chrY:
            chromosomes.discard("chrY")
    #            chromosomes.discard('chr9')
    else:
        chromosomes = args.chromosomes.split(",")

    chrom_sizes_map = pd.read_csv(
        args.chrom_sizes, sep="\t", header=None, index_col=0
    ).to_dict()[1]

    powerlaw_params = pd.read_csv(args.powerlaw_params_tsv, sep="\t").iloc[0]
    hic_gamma, hic_scale = powerlaw_params["hic_gamma"], powerlaw_params["hic_scale"]
    for chromosome in chromosomes:
        print("Making predictions for chromosome: {}".format(chromosome))
        t = time.time()
        this_enh = enhancers.loc[enhancers["chr"] == chromosome, :].copy()
        this_genes = genes.loc[genes["chr"] == chromosome, :].copy()

        this_chr = make_predictions(
            chromosome,
            this_enh,
            this_genes,
            args,
            hic_gamma,
            hic_scale,
            chrom_sizes_map,
        )
        all_putative_list.append(this_chr)

        print(
            "Completed chromosome: {}. Elapsed time: {} \n".format(
                chromosome, time.time() - t
            )
        )

    # Subset predictions
    print("Writing output files...")
    all_putative = pd.concat(all_putative_list)
    all_putative["CellType"] = args.cellType
    if args.hic_dir:
        all_putative["hic_contact_squared"] = all_putative["hic_contact"] ** 2
    slim_cols = [
        "chr",
        "start",
        "end",
        "name",
        "TargetGene",
        "TargetGeneTSS",
        "CellType",
        args.score_column,
    ]
    if args.run_all_genes:
        all_positive = all_putative.iloc[
            np.logical_and.reduce(
                (
                    all_putative[args.score_column] > args.threshold,
                    ~(all_putative["class"] == "promoter"),
                )
            ),
            :,
        ]
    else:
        all_positive = all_putative.iloc[
            np.logical_and.reduce(
                (
                    all_putative.TargetGeneIsExpressed,
                    all_putative[args.score_column] > args.threshold,
                    ~(all_putative["class"] == "promoter"),
                )
            ),
            :,
        ]

    if args.include_self_promoter:
        self_promoter = all_putative.loc[all_putative["isSelfPromoter"], :]
        all_positive = pd.concat([all_positive, self_promoter]).drop_duplicates()

    all_positive.to_csv(
        pred_file_full, sep="\t", index=False, header=True, float_format="%.6f"
    )
    all_positive[slim_cols].to_csv(
        pred_file_slim, sep="\t", index=False, header=True, float_format="%.6f"
    )

    make_gene_prediction_stats(all_putative, args)
    write_connections_bedpe_format(all_positive, pred_file_bedpe, args.score_column)

    if args.make_all_putative:
        if not args.use_hdf5:
            all_putative.loc[all_putative.TargetGeneIsExpressed, :].to_csv(
                all_pred_file_expressed,
                sep="\t",
                index=False,
                header=True,
                compression="gzip",
                float_format="%.6f",
                na_rep="NaN",
            )
            all_putative.loc[~all_putative.TargetGeneIsExpressed, :].to_csv(
                all_pred_file_nonexpressed,
                sep="\t",
                index=False,
                header=True,
                compression="gzip",
                float_format="%.6f",
                na_rep="NaN",
            )
        else:
            all_pred_file_expressed = os.path.join(
                args.outdir, "EnhancerPredictionsAllPutative.h5"
            )
            all_pred_file_nonexpressed = os.path.join(
                args.outdir, "EnhancerPredictionsAllPutativeNonExpressedGenes.h5"
            )
            all_putative.loc[all_putative.TargetGeneIsExpressed, :].to_hdf(
                all_pred_file_expressed, key="predictions", complevel=9, mode="w"
            )
            all_putative.loc[~all_putative.TargetGeneIsExpressed, :].to_hdf(
                all_pred_file_nonexpressed, key="predictions", complevel=9, mode="w"
            )

    test_variant_overlap(args, all_putative)

    print("Done.")


def validate_args(args):
    if args.hic_dir and args.hic_type == "juicebox":
        assert (
            args.hic_resolution is not None
        ), "HiC resolution must be provided if hic_type is juicebox"

    if not args.hic_dir:
        print(
            "WARNING: Hi-C directory not provided. Model will only compute ABC score using powerlaw!"
        )


if __name__ == "__main__":
    main()
