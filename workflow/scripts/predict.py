import argparse
from predictor import *
from tools import *
from getVariantOverlap import *
import pandas as pd
import numpy as np
import sys, traceback, os, os.path
import time

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
    parser.add_argument('--score_column', default='ABC.Score', help="Column name of score to use for thresholding")
    parser.add_argument('--threshold', type=float, required=True, default=.022, help="Threshold on ABC Score (--score_column) to call a predicted positive")
    parser.add_argument('--cellType', help="Name of cell type")
    parser.add_argument('--chrom_sizes', help="Chromosome sizes file")

    #hic
    #To do: validate params
    parser.add_argument('--HiCdir', default=None, help="HiC directory")
    parser.add_argument('--hic_resolution', type=int, help="HiC resolution")
    parser.add_argument('--tss_hic_contribution', type=float, default=100, help="Weighting of diagonal bin of hic matrix as a percentage of the maximum of its neighboring bins")
    parser.add_argument('--hic_pseudocount_distance', type=int, default=1e6, help="A pseudocount is added equal to the powerlaw fit at this distance")
    parser.add_argument('--hic_type', default = 'juicebox', choices=['juicebox','bedpe', 'avg'], help="format of hic files")
    parser.add_argument('--hic_is_doubly_stochastic', action='store_true', help="If hic matrix is already DS, can skip this step")

    #Power law
    parser.add_argument('--scale_hic_using_powerlaw', action="store_true", help="Scale Hi-C values using powerlaw relationship")
    parser.add_argument('--hic_gamma', type=float, default=.87, help="powerlaw exponent of hic data. Must be positive")
    parser.add_argument('--hic_scale', type=float, help="scale of hic data. Must be positive")
    parser.add_argument('--hic_gamma_reference', type=float, default=.87, help="powerlaw exponent to scale to. Must be positive")

    #Genes to run through model
    parser.add_argument('--run_all_genes', action='store_true', help="Do not check for gene expression, make predictions for all genes")
    parser.add_argument('--expression_cutoff', type=float, default=1, help="Make predictions for genes with expression higher than this value")
    parser.add_argument('--promoter_activity_quantile_cutoff', type=float, default=.4, help="Quantile cutoff on promoter activity. Used to consider a gene 'expressed' in the absence of expression data")

    #Output formatting
    parser.add_argument('--make_all_putative', action="store_true", help="Make big file with concatenation of all genes file")
    parser.add_argument('--use_hdf5', action="store_true", help="Write AllPutative file in hdf5 format instead of tab-delimited")

    #Other
    parser.add_argument('--tss_slop', type=int, default=500, help="Distance from tss to search for self-promoters")
    parser.add_argument('--chromosomes', default="all", help="chromosomes to make predictions for. Defaults to intersection of all chromosomes in --genes and --enhancers")
    parser.add_argument('--include_chrY', '-y', action='store_true', help="Make predictions on Y chromosome")

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
    genes = pd.read_csv(args.genes, sep = "\t")
    genes = determine_expressed_genes(genes, args.expression_cutoff, args.promoter_activity_quantile_cutoff)
    genes = genes.loc[:,['chr','symbol','tss','Expression','PromoterActivityQuantile','isExpressed', 'ATAC.RPKM.quantile.TSS1Kb']]
    genes.columns = ['chr','TargetGene', 'TargetGeneTSS', 'TargetGeneExpression', 'TargetGenePromoterActivityQuantile','TargetGeneIsExpressed', 'normalized_atac']

    print("reading enhancers")
    enhancers_full = pd.read_csv(args.enhancers, sep = "\t")
    #TO DO
    #Think about which columns to include
    enhancers = enhancers_full.loc[:,['chr','start','end','name','class','activity_base','normalized_atac' ]]
    enhancers['activity_base_squared'] = enhancers['activity_base']**2
    #Initialize Prediction files
    pred_file_full = os.path.join(args.outdir, "EnhancerPredictionsFull.txt")
    pred_file_slim = os.path.join(args.outdir, "EnhancerPredictions.txt")
    pred_file_bedpe = os.path.join(args.outdir, "EnhancerPredictions.bedpe")
    all_pred_file_expressed = os.path.join(args.outdir, "EnhancerPredictionsAllPutative.txt.gz")
    all_pred_file_nonexpressed = os.path.join(args.outdir, "EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz")
    variant_overlap_file = os.path.join(args.outdir, "EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.txt.gz")
    all_putative_list = []

    #Make predictions
    if args.chromosomes == "all":
        chromosomes = set(genes['chr']).intersection(set(enhancers['chr'])) 
        if not args.include_chrY:
            chromosomes.discard('chrY')
#            chromosomes.discard('chr9')
    else:
        chromosomes = args.chromosomes.split(",")

    for chromosome in chromosomes:
        print('Making predictions for chromosome: {}'.format(chromosome))
        t = time.time()
        this_enh = enhancers.loc[enhancers['chr'] == chromosome, :].copy()
        this_genes = genes.loc[genes['chr'] == chromosome, :].copy()

        this_chr = make_predictions(chromosome, this_enh, this_genes, args)
        all_putative_list.append(this_chr)

        print('Completed chromosome: {}. Elapsed time: {} \n'.format(chromosome, time.time() - t))

    # Subset predictions
    print("Writing output files...")
    all_putative = pd.concat(all_putative_list)
    all_putative['CellType'] = args.cellType
    all_putative['hic_contact_squared'] = all_putative['hic_contact']**2
    slim_cols = ['chr','start','end','name','TargetGene','TargetGeneTSS','CellType',args.score_column]
    if args.run_all_genes:
        all_positive = all_putative.iloc[np.logical_and.reduce((all_putative[args.score_column] > args.threshold, ~(all_putative['class'] == "promoter"))),:]
    else:
        all_positive = all_putative.iloc[np.logical_and.reduce((all_putative.TargetGeneIsExpressed, all_putative[args.score_column] > args.threshold, ~(all_putative['class'] == "promoter"))),:]

    all_positive.to_csv(pred_file_full, sep="\t", index=False, header=True, float_format="%.6f")
    all_positive[slim_cols].to_csv(pred_file_slim, sep="\t", index=False, header=True, float_format="%.6f")

    
    make_gene_prediction_stats(all_putative, args)
    write_connections_bedpe_format(all_positive, pred_file_bedpe, args.score_column)
    
    if args.make_all_putative:
        if not args.use_hdf5:
            all_putative.loc[all_putative.TargetGeneIsExpressed,:].to_csv(all_pred_file_expressed, sep="\t", index=False, header=True, compression="gzip", float_format="%.6f", na_rep="NaN")
            all_putative.loc[~all_putative.TargetGeneIsExpressed,:].to_csv(all_pred_file_nonexpressed, sep="\t", index=False, header=True, compression="gzip", float_format="%.6f", na_rep="NaN")
        else:
            all_pred_file_expressed = os.path.join(args.outdir, "EnhancerPredictionsAllPutative.h5")
            all_pred_file_nonexpressed = os.path.join(args.outdir, "EnhancerPredictionsAllPutativeNonExpressedGenes.h5")
            all_putative.loc[all_putative.TargetGeneIsExpressed,:].to_hdf(all_pred_file_expressed, key='predictions', complevel=9, mode='w')
            all_putative.loc[~all_putative.TargetGeneIsExpressed,:].to_hdf(all_pred_file_nonexpressed, key='predictions', complevel=9, mode='w')
   
    test_variant_overlap(args, all_putative.loc[all_putative.TargetGeneIsExpressed,:])
    print("Done.")
    
def validate_args(args):
    if args.HiCdir and args.hic_type == 'juicebox':
        assert args.hic_resolution is not None, 'HiC resolution must be provided if hic_type is juicebox'

    if not args.HiCdir:
        print("WARNING: Hi-C directory not provided. Model will only compute ABC score using powerlaw!")

if __name__ == '__main__':
    main()
    
