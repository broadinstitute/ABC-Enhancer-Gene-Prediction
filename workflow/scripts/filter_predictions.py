import click
import pandas as pd
from predictor import make_gene_prediction_stats
from tools import write_connections_bedpe_format


@click.command()
@click.option("--output_tsv_file", type=str)
@click.option("--output_slim_tsv_file", type=str)
@click.option("--output_bed_file", type=str)
@click.option("--output_gene_stats_file", type=str)
@click.option("--pred_file", type=str)
@click.option("--pred_nonexpressed_file", type=str)
@click.option(
    "--score_column",
    type=str,
    help="Column name of score to use for thresholding",
)
@click.option("--threshold", type=float)
@click.option("--include_self_promoter", type=bool)
@click.option("--only_expressed_genes", type=bool)
def main(
    output_tsv_file,
    output_slim_tsv_file,
    output_bed_file,
    output_gene_stats_file,
    pred_file,
    pred_nonexpressed_file,
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

    write_connections_bedpe_format(filtered_predictions, output_bed_file, score_column)
    make_gene_prediction_stats(
        filtered_predictions, score_column, threshold, output_gene_stats_file
    )


def remove_promoters(pred_df: pd.DataFrame, keep_self_promoters: bool) -> pd.DataFrame:
    if keep_self_promoters:
        return pred_df[(pred_df["class"] != "promoter") | pred_df["isSelfPromoter"]]
    else:
        return pred_df[pred_df["class"] != "promoter"]


if __name__ == "__main__":
    main()
