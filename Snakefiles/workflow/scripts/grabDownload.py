#! bin/bash 

import sys, os
import argparse 
import subprocess
import numpy as np
import pandas as pd
from utils import *
from subprocess import check_call, check_output, PIPE, Popen, getoutput, CalledProcessError

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dhs', help="DNase experiments file")
    parser.add_argument('--h3k27ac', help="H3K27ac experiments file")
    parser.add_argument('--atac', help="ATAC experiments file")
    parser.add_argument('--genome_assembly', help="Genome assembly")
    parser.add_argument('--download_files', action='store_true', help="Flag to download files or not")
    parser.add_argument('--outdir', default=".", help="Outdir to save downloaded metadata files")
    parser.add_argument('--data_outdir', help="Outdir to save downloaded bam files")
    args = parser.parse_args()
    return args 

def main():
    args = parse_args()
    metadata = assignFiltersToDataFrame(args)
    if args.download_files:
        outfile = downloadFiles(args, metadata)
    save_paired_single_end_files(args, metadata)


def assignFiltersToDataFrame(args):
   
    # TODO: Filters work on current hg19 file
    # Add filters to deal with paired end, single end
    # Paired end flag in metadata doesn't seem to be informative 

    # open dataframes 
    rint("Opening DHS and H3K27ac Experiment Files...")
    dhs_data = pd.read_csv(args.dhs, sep="\t")
    h3k27ac_data = pd.read_csv(args.h3k27ac, sep="\t")

    print("Filtering for genome assembly...")
    # filter for assembly hg19
    # fill up dhs mapped read length table with numbers
    

    dhs_alignment_bam = mapExperimentToLength(dhs_data)
    # H3K27ac goes through a different filtering because 
    h3k27ac_alignment_bam = mapExperimentToLength(h3k27ac_data)

    merge_columns = ['Biosample term name','Biosample organism',  'Biosample treatments','Biosample treatments amount', 'Biosample treatments duration','Biosample genetic modifications methods','Biosample genetic modifications categories','Biosample genetic modifications targets','Biosample genetic modifications gene targets', 'Assembly', 'Biosample type', 'File Status', 'File format', 'File type', 'Library made from', 'Output type', 'Biosample term id', 'Genome annotation', 'File format type', 'Output type', 'Mapped read length']
   
    # merge dhs and h3k27ac
    intersected = pd.merge(dhs_alignment_bam, h3k27ac_alignment_bam, how='inner', on=merge_columns, suffixes=('_Accessibility', '_H3K27ac'))
    
    if args.atac is not None:
        atac_df = pd.read_csv(args.atac, sep="\t")
        atac = mapExperimentToLength(atac_df)
        intersected_df = pd.merge(atac, h3k27ac_data, how='inner', on=merge_columns, suffixes=('_Accessibility', '_H3K27ac'))
        copy = intersected
        intersected = pd.concat([copy, intersected_df])
    # filter for filtered file + released files 
    #filtered_intersected = intersected.loc[intersected['Output type'].eq('alignments') & intersected['File Status'].eq('released')]
    # fill columns that are filled with NAN
    df = intersected.fillna(0.0)
    
    # grab entries with biological replicates 
    # grab celltypes with biological replicates 
    full_metadata, metadata_unique = obtainDuplicated(args, df)
    # grab entries with no biological repliates 
    # filter for mapped read lengths of usually 32.0 or 36.0

    # metadata_unique contains the combined accession numbers for experiments with technical/biological replicates 
    # full metadata contains every technical and biological replicate as its own single entry
    # save relevant columns into input data lookup for input into ABC code
    prepareLookup(args, metadata_unique, "input_data_lookup")
    return full_metadata

def downloadFiles(args, df):
    # get download links and start downloading 
    download_links = df[['File download URL_Accessibility', 'File download URL_H3K27ac']].melt(value_name='download_links')
    outfile=os.path.join(args.outdir, "linkstodownload.txt")
    download_links[['download_links']].to_csv( outfile, sep="\t", index=False, header=None)
    
    if not os.path.exists(args.data_outdir):
        os.mkdir(args.data_outdir)

    commands = "bash downloadFiles.sh {} {}".format(outfile, args.data_outdir)
    process = subprocess.call(commands, shell=True)
    return outfile 

def save_paired_single_end_files(args, metadata):
    for col, title in zip(["single-ended", "paired-ended"], ["singleend", "pairedend"]):
        paired = metadata['File accession_Accessibility'].loc[metadata['Run type_Accessibility']==str(col)].drop_duplicates()
        paired_h3k27ac = metadata['File accession_H3K27ac'].loc[metadata['Run type_H3K27ac']==str(col)].drop_duplicates()
        combined_paired = pd.concat([paired, paired_h3k27ac])
        combined_paired.to_csv(os.path.join(args.outdir, "unique_{}_h3k27ac_dhs_files.tsv".format(str(title))), sep="\t", index=False)


if __name__=="__main__":
    main()
