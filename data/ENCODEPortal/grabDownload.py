#! bin/bash 

import sys, os
import argparse 
import subprocess
import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dhs', help="DNase experiments file")
    parser.add_argument('--h3k27ac', help="H3K27ac experiments file")
    parser.add_argument('--atac', help="ATAC experiments file")
    parser.add_argument('--genome_assembly', help="Genome assembly")
    parser.add_argument('--include_extra_map_length', default=False, help="Includes entries with Mapped Read Lengths > 40")
    parser.add_argument('--download_files', default=False, help="Flag to download files or not")
    parser.add_argument('--outdir', default=".", help="Outdir to save downloaded bam files")
    args = parser.parse_args()
    return args 

def main():
    args = parse_args()
    df = assignFiltersToDataFrame(args)
    if args.download_files:
        outfile = downloadFiles(df, args.outdir)

def assignFiltersToDataFrame(args):
    
    # TODO: Filters work on current hg19 file
    # Add filters to deal with paired end, single end

    # open dataframes 
    print("Opening DHS and H3K27ac Experiment Files...")
    dhs_data = pd.read_csv(args.dhs, sep="\t")
    h3k27ac_data = pd.read_csv(args.h3k27ac, sep="\t")

    print("Filtering for genome assembly...")
    # filter for assembly hg19
    dhs = dhs_data.loc[dhs_data['Assembly'] == args.genome_assembly]
    h3k27ac = h3k27ac_data.loc[h3k27ac_data['Assembly'] == args.genome_assembly]

    merge_columns = ['Biosample term name','Biosample organism', 'Biosample treatments','Biosample treatments amount', 'Biosample treatments duration','Biosample genetic modifications methods','Biosample genetic modifications categories','Biosample genetic modifications targets','Biosample genetic modifications gene targets']
    
    # merge dhs and h3k27ac
    intersected = pd.merge(dhs, h3k27ac, how='inner', on=merge_columns)
    
    if args.atac is not None:
        atac = pd.read_csv(args.atac, sep="\t")
        intersected_df = intersected
        intersected = pd.merge(intersected_df, atac, how='inner', on=merge_columns)

    # filter for filtered file + released files 
    filtered_intersected = intersected.loc[intersected['Output type_x'].eq('alignments') & intersected['Output type_y'].eq('alignments') & intersected['File Status_x'].eq('released') & intersected['File Status_y'].eq('released')]
    # remove columns that are filled with NAN
    df = filtered_intersected.fillna(0.0)
    subset_intersected = df.loc[:, (df != 0).all(axis=0)]

    # grab entries with biological replicates 
    duplicates = subset_intersected[subset_intersected.duplicated(['Experiment accession_x'], keep=False) & ~subset_intersected.duplicated(['Experiment accession_x', 'Biological replicate(s)_x'], keep=False)].drop_duplicates(['Biosample term name', 'Biological replicate(s)_x'])
    # grab celltypes with biological replicates 
    df_biological_rep = duplicates['Biosample term name'].drop_duplicates()

    # grab entries with no biological replicates 
    # filter for mapped read lengths of usually 32.0 or 36.0
    columns = ['Biosample term name']
    final_experiment_file = subset_intersected.loc[np.logical_not(subset_intersected['Biosample term name'].isin(df_biological_rep)) & subset_intersected['Mapped read length_y'].between(30.0,40.0)].drop_duplicates(columns)
    
    combined_experiments = pd.concat([final_experiment_file, duplicates])
    # include entries with higher mapped read lengths
    if args.include_extra_map_length:
        mapped_experiment_file = subset_intersected.loc[np.logical_not(subset_intersected['Biosample term name'].isin(df_biological_rep)) & subset_intersected['Mapped read length_y'].between(40.0, 200.0)].drop_duplicates(columns)
        df_2 = combined_experiments
        combined_experiments = pd.concat([df_2, mapped_experiment_file])

    combined_experiments.to_csv("metadata.tsv", sep="\t")
    return combined_experiments

def downloadFiles(df, outdir):
    # get download links and start downloading 
    download_links = df[['File download URL_x', 'File download URL_y']].melt(value_name='download_links')
    outfile = "linkstodownload.txt"
    download_links[['download_links']].to_csv(outfile, sep="\t", index=False, header=None)
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    commands = "parallel --verbose -j 10 bash downloadFiles.sh ::: {} :::  {}".format(outfile, outdir)
    process = subprocess.call(commands, shell=True)
    return outfile 


if __name__=="__main__":
    main()
