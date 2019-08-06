#!/usr/bin/env python3

import argparse
import progressbar as pb
from predictor import Predictor
from tools import *
import pandas as pd
import numpy as np
import sys, traceback, os, os.path

def get_model_argument_parser():
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description='Predict enhancer relative effects.',
                                     formatter_class=formatter)
    readable = argparse.FileType('r')

    #Basic parameters
    parser.add_argument('--enhancers', required=True, help="Candidate enhancer regions. Formatted as the EnhancerList.txt file produced by run.neighborhoods.py")
    parser.add_argument('--genes', required=True, help="Genes to make predictions for. Formatted as the GeneList.txt file produced by run.neighborhoods.py")
    parser.add_argument('--outdir', required=True, help="output directory")
    parser.add_argument('--window', type=int, default=5000000, help="Make predictions for all candidate elements within this distance of the gene's TSS")
    parser.add_argument('--threshold', type=float, required=True, default=None, help="Threshold on ABC Score to call a predicted positive")
    parser.add_argument('--cellType', help="Name of cell type")

    #hic
    parser.add_argument('--HiCdir', help="Directory with hic bedgraphs")
    parser.add_argument('--tss_hic_contribution', type=float, default=100, help="Weighting of diagonal bin of hic matrix as a percentage of its neighbors")
    parser.add_argument('--hic_pseudocount_distance', type=int, default=1e6, help="A pseudocount is added equal to the powerlaw fit at this distance")

    #Power law
    parser.add_argument('--scale_hic_using_powerlaw', action="store_true", help="Scale Hi-C values using powerlaw relationship")
    parser.add_argument('--hic_gamma', type=float, default=1, help="powerlaw exponent of hic data. Must be positive")
    parser.add_argument('--hic_gamma_reference', type=float, default=1, help="powerlaw exponent to scale to. Must be positive")

    #Genes to run through model
    parser.add_argument('--run_all_genes', action='store_true', help="Do not check for gene expression, make predictions for all genes")
    parser.add_argument('--expression_cutoff', type=float, default=1, help="Make predictions for genes with expression higher than this value")
    parser.add_argument('--promoter_activity_quantile_cutoff', type=float, default=.4, help="Quantile cutoff on promoter activity. Used to consider a gene 'expressed' in the absence of expression data")

    #Output formatting
    parser.add_argument('--skip_gene_files', action="store_true", help="Do not make individual gene files")
    parser.add_argument('--skinny_gene_files', action="store_true", help="Use subset of columns for genes files")
    parser.add_argument('--make_all_putative', action="store_true", help="Make big file with concatenation of all genes file")

    #Other
    parser.add_argument('--tss_slop', type=int, default=500, help="Distance from tss to search for self-promoters")
    parser.add_argument('--include_chrY', '-y', action='store_true', help="Include Y chromosome")

    return parser


def get_predict_argument_parser():
    parser = get_model_argument_parser()
    return parser

def main():
    parser = get_predict_argument_parser()
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    preddir = os.path.join(args.outdir, "genes")
    if not os.path.exists(preddir):
        os.makedirs(preddir)

    write_prediction_params(args, os.path.join(args.outdir, "parameters.predict.txt"))
    
    print("reading genes")
    genes = read_genes(args.genes)

    print("reading enhancers")
    enhancers = read_enhancers(args.enhancers)

    print("building predictor")
    predictor = Predictor(enhancers, **vars(args))

    chromosomes = predictor.chromosomes()
    print("data loaded for chromosomes: {}".format(" ".join(sorted(chromosomes))))

    #Initialize Prediction files
    pred_file = os.path.join(args.outdir, "EnhancerPredictions.txt")
    all_pred_file = os.path.join(args.outdir, "EnhancerPredictionsAllPutative.txt.gz")

    args.score_column = "ABC.Score"

    all_positive_list = []
    all_putative_list = []
    gene_stats = []
    failed_genes = []

    for idx, gene in genes.iterrows():
        if gene.chr == 'chrY' and not args.include_chrY:
            continue
        if gene.chr not in chromosomes:
            print("\nNo data for {}".format(gene.chr))
            continue
        print("\nPredicting {} with TSS {} {}".format(gene["name"], gene["chr"], gene["tss"]))

        try:
            nearby_enhancers = enhancers.within_range(gene.chr, gene.tss - args.window, gene.tss + args.window)
            predictor.predict_from_normalized_to_enhancers(nearby_enhancers, gene, args.window, tss_slop=args.tss_slop)
            
            col_names = ['chr','start','end','class','TargetGene','TargetGeneTSS','ABC.Score','powerlaw.Score','distance','hic.distance','hic.distance.adj','estimatedCP','estimatedCP.adj','normalized_dhs','normalized_atac','activity_base','normalized_h3K27ac','TargetGeneExpression','TargetGenePromoterActivityQuantile']
            col_names = [col for col in col_names if col in nearby_enhancers.columns]
            #col_names = list(set(nearby_enhancers.columns) & set(col_names))

            if not args.skip_gene_files:
                if not args.skinny_gene_files:
                    write_scores(preddir, gene, nearby_enhancers)
                else:
                    write_scores(preddir, gene, nearby_enhancers[col_names])

            if args.make_all_putative:
                if args.skinny_gene_files:
                    all_putative_list.append(nearby_enhancers[col_names])
                else:
                    all_putative_list.append(nearby_enhancers)
            
            gene_is_expressed_proxy = check_gene_for_runnability(gene, args.expression_cutoff, args.promoter_activity_quantile_cutoff)
            if args.run_all_genes or gene_is_expressed_proxy:
                positives = nearby_enhancers.ix[np.logical_and(nearby_enhancers[args.score_column] >= args.threshold, nearby_enhancers["class"] != "promoter"),:]
                print("{} enhancers predicted for {}".format(positives.shape[0], gene["name"]))
                all_positive_list.append(positives)

            #Consider a gene as failed if all its ABC Scores are nan
            if all(nearby_enhancers[args.score_column].isnull().tolist()):
                failed_genes.append(gene["chr"] + "\t" + gene["name"])

            #Add gene to gene summary file
            if nearby_enhancers.shape[0] > 0:
                stats = predictor.get_gene_prediction_stats(args, nearby_enhancers)
                stats['prediction_file'] = get_score_filename(gene)
                stats['gene_is_expressed_proxy'] = gene_is_expressed_proxy
                gene_stats.append(stats)

        except:
            failed_genes.append(gene["chr"] + "\t" + gene["name"])
            print("Failed on " + gene["name"] + " ... skipping. Traceback:")
            traceback.print_exc(file=sys.stdout)
            continue

    if args.score_column is not None:
        all_positive = pd.concat(all_positive_list)
        all_positive.to_csv(pred_file, sep="\t", index=False, header=True, float_format="%.4f")
        write_connections_bedpe_format(all_positive, outfile=os.path.join(args.outdir, "EnhancerPredictions.bedpe"), score_column=args.score_column)

    gene_stats = pd.concat(gene_stats, axis=1).T
    gene_stats.to_csv(os.path.join(args.outdir, "GenePredictionStats.txt"), sep="\t", index=False)

    with open(os.path.join(args.outdir, "FailedGenes.txt"), 'w') as failed_file:
        for gene in failed_genes:
            failed_file.write(gene + "\n")

    if args.make_all_putative:
        all_putative = pd.concat(all_putative_list)
        all_putative.to_csv(all_pred_file, sep="\t", index=False, header=True, compression="gzip", float_format="%.4f", na_rep="NaN", chunksize=100000)

def write_prediction_params(args, file):
    with open(file, 'w') as outfile:
        for arg in vars(args):
            outfile.write("--" + arg + " " + str(getattr(args, arg)) + " ")


if __name__ == '__main__':
    main()
