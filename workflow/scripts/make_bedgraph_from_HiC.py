import argparse
import glob
import os.path
from hic import HiC
import pandas
import numpy as np
import sys
from neighborhoods import read_bed, process_gene_bed

def parseargs():
    parser = argparse.ArgumentParser(description='Convert HiC matrices to bedgraphs for a set of genes')
    parser.add_argument('--outdir', required=True, help="directory to write HiC bedgraphs")
    parser.add_argument('--hic_dir', required=True, help="location of HiC data")
    
    #Genes    
    parser.add_argument('--genes', required=True, help=".bed file of genes")
    parser.add_argument('--gene_name_annotations', default="symbol", help="Comma delimited string of names corresponding to the gene identifiers present in the name field of the gene annotation bed file")
    parser.add_argument('--primary_gene_identifier', default="symbol", help="Primary identifier used to identify genes. Must be present in gene_name_annotations")

    #HiC Params
    parser.add_argument('--resolution', type=int, default=5000, help="HiC resolution to use")
    parser.add_argument('--kr_cutoff', type=float, default=0.1, help="Measured data from Hi-C matrix for rows/columns with kr normalization vector below this value are not used. Instead they are interpolated from neighboring bins")
    parser.add_argument('--window', type=int, default=5000000, help="maximum distance from each TSS to store (bp)")

    parser.add_argument('--overwrite', action="store_true", help="force overwriting files")

    return parser.parse_args()


if __name__ == '__main__':
    args = parseargs()

    def match_files(*subdirs):
        return glob.glob(os.path.join(args.hic_dir, *subdirs))

    #Read genes
    genes_bed = read_bed(args.genes) 
    genes = process_gene_bed(genes_bed, args.gene_name_annotations, args.primary_gene_identifier, fail_on_nonunique=False)

    #Get raw hic and normalization files
    resolution = '{}kb'.format(args.resolution // 1000)
    hic_files = {}
    for chr in set(genes['chr']):
        possible_files = match_files('{}'.format(chr), '{}_{}.RAWobserved'.format(chr, resolution))
        possible_norms = match_files('{}'.format(chr), '{}_{}.KRnorm'.format(chr, resolution))

        if possible_files:
            hic_files['{}'.format(chr)] = (possible_files[0], possible_norms[0])

    # create data accessor
    hic_data = HiC(hic_files, window=args.window, resolution=args.resolution, kr_cutoff=args.kr_cutoff)

    # create output directory
    os.makedirs(args.outdir, exist_ok=True)
    #os.makedirs(os.path.join(args.outdir, "raw"), exist_ok=True)

    #Make a bedgraph per gene
    skipped = []
    for idx, gene in genes.iterrows():
        if gene.chr not in hic_data.chromosomes():
            print("No HiC data for {} on {}".format(gene['name'], gene.chr))
            continue
        filename = os.path.join(args.outdir,
                        "{}_{}_{}.bg.gz".format(gene['name'] or "UNK", gene.chr, int(gene.tss)))

        if not args.overwrite:
            if os.path.exists(filename):
                skipped.append(gene['name'])
                print("Skipping {} on {} with tss {} since it already has hic data and --overwrite flag is not set".format(gene['name'], gene.chr, gene.tss))
                continue

        hic_row = hic_data.row(gene.chr, gene.tss)

        #Include all values within args.window of the tss. This will facilitate interpolating NaNs
        values = [(gene.chr, idx * args.resolution,
               (idx + 1) * args.resolution,
               hic_row[0, idx]) for idx in range(hic_row.A.shape[1]) if abs(idx * args.resolution - int(gene.tss)) < args.window]

        #interpolate the nan's. Sometimes there may be nan's at the beginning/end of the vector - set to 0
        #note this is only interpreting the nan's: missing data due to low kr norm value. This is not interpolating 0's in HiC
        df2 = pandas.DataFrame.from_records(values).interpolate().fillna(value=0)
        df2.to_csv(filename,
              sep='\t', compression='gzip',
              header=False, index=False)

        print("Completed {} on {}".format(gene['name'], gene.chr))

    if len(skipped) > 0:
        print("Skipped {} genes because they already have HiC files".format(len(skipped)))
