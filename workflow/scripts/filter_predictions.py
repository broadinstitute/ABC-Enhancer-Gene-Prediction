import click
import argparse
import pandas as pd
from predictor import make_gene_prediction_stats
from tools import write_connections_bedpe_format
from getVariantOverlap import test_variant_overlap


@click.command()
@click.option("--outdir", required=True, help="output directory")
@click.option("--output_tsv_file", type=str)
@click.option("--output_slim_tsv_file", type=str)
@click.option("--output_bed_file", type=str)
@click.option("--output_gene_stats_file", type=str)
@click.option("--pred_file", type=str)
@click.option("--pred_nonexpressed_file", type=str)
@click.option("--chrom_sizes", type=str)
@click.option("--score_column", type=str,help="Column name of score to use for thresholding")
@click.option("--threshold", type=float)
@click.option("--include_self_promoter", type=bool)
@click.option("--only_expressed_genes", type=bool)

def main(
    outdir,
    output_tsv_file,
    output_slim_tsv_file,
    output_bed_file,
    output_gene_stats_file,
    pred_file,
    pred_nonexpressed_file,
    chrom_sizes,
    score_column,
    threshold,
    include_self_promoter,
    only_expressed_genes,
):
    slim_columns = [
        "chr",
        "start",
        "end",
        "name",
        "TargetGene",
        "TargetGeneTSS",
        "CellType",
        score_column,
    ]
    print("debug1")
    all_putative = pd.read_csv(pred_file, sep="\t")
    if not only_expressed_genes:
        non_expressed = pd.read_csv(pred_nonexpressed_file, sep="\t")
        all_putative = pd.concat([all_putative, non_expressed], ignore_index=True)

    filtered_predictions = all_putative[all_putative[score_column] > threshold]
    filtered_predictions = remove_promoters(filtered_predictions, include_self_promoter)

    filtered_predictions.to_csv(
        output_tsv_file, sep="\t", index=False, header=True, float_format="%.6f"
    )
    filtered_predictions_slim = filtered_predictions[slim_columns]
    filtered_predictions_slim.to_csv(
        output_slim_tsv_file, sep="\t", index=False, header=True, float_format="%.6f"
    )
    print("debug2")
    parser=get_variant_overlapping_argument_parser()
    print(parser)
    args=parser.parse_args([
    "--score_column", score_column,
    "--chrom_sizes", chrom_sizes,
    "--outdir", outdir,
    "--threshold", str(threshold)
    ])
    print(args)
    test_variant_overlap(args, all_putative)

    write_connections_bedpe_format(filtered_predictions, output_bed_file, score_column)
    make_gene_prediction_stats(
        filtered_predictions, score_column, threshold, output_gene_stats_file
    )


def get_model_argument_parser():
    class formatter(
        argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
    ):
        pass
    parser = argparse.ArgumentParser(
        description="Filter and shrink putative enhancers for variant overlapping.", formatter_class=formatter
    )
    readable = argparse.FileType("r")
    parser.add_argument("--score_column",default="ABC.Score",help="Column name of score to use for thresholding")
    parser.add_argument("--chrom_sizes", required=True, help="Chromosome sizes file")
    parser.add_argument("--outdir", required=True, help="output directory")
    parser.add_argument("--threshold", required=True, help="threshold")
    return parser

def get_variant_overlapping_argument_parser():
    parser = get_model_argument_parser()
    return parser

def remove_promoters(pred_df: pd.DataFrame, keep_self_promoters: bool) -> pd.DataFrame:
    if keep_self_promoters:
        return pred_df[(pred_df["class"] != "promoter") | pred_df["isSelfPromoter"]]
    else:
        return pred_df[pred_df["class"] != "promoter"]


if __name__ == "__main__":
    main()
