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
    parser.add_argument('--include_extra_map_length', action='store_true', help="Includes entries with Mapped Read Lengths > 40")
    parser.add_argument('--download_files', action='store_true', help="Flag to download files or not")
    parser.add_argument('--outdir', default=".", help="Outdir to save downloaded bam files")
    args = parser.parse_args()
    return args 

def main():
    args = parse_args()
    print(args)
    metadata = assignFiltersToDataFrame(args)
    if args.download_files:
        print(args.download_files)
        outfile = downloadFiles(args, metadata)
    get_paired_single_end_datafiles(args, metadata)


def assignFiltersToDataFrame(args):
   
    # TODO: Filters work on current hg19 file
    # Add filters to deal with paired end, single end
    # Paired end flag in metadata doesn't seem to be informative 

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
        intersected_df = pd.merge(atac, h3k27ac, how='inner', on=merge_columns)
        copy = intersected
        intersected = pd.concat([copy, intersected_df])

    # filter for filtered file + released files 
    filtered_intersected = intersected.loc[intersected['Output type_x'].eq('alignments') & intersected['Output type_y'].eq('alignments') & intersected['File Status_x'].eq('released') & intersected['File Status_y'].eq('released')]
    # remove columns that are filled with NAN
    df = filtered_intersected.fillna(0.0)
    subset_intersected = df.loc[:, (df != 0).all(axis=0)]

    # grab entries with biological replicates 
    # grab celltypes with biological replicates 
    duplicates, df_biological_rep, metadata_unique = obtainDuplicated(args, subset_intersected)
    # grab entries with no biological replicates 
    # filter for mapped read lengths of usually 32.0 or 36.0
    columns = ['Biosample term name']
    final_experiment_file = subset_intersected
    final_experiment_file = subset_intersected.loc[np.logical_not(subset_intersected['Biosample term name'].isin(df_biological_rep))].drop_duplicates(columns)
#    final_experiment_file = subset_intersected.loc[np.logical_not(subset_intersected['Biosample term name'].isin(df_biological_rep)) & subset_intersected['Mapped read length_y'].between(30.0,40.0)].drop_duplicates(columns)
    combined_experiments = pd.concat([final_experiment_file, duplicates])
    # include entries with higher mapped read lengths
#    if args.include_extra_map_length:
#        mapped_experiment_file = subset_intersected.loc[np.logical_not(subset_intersected['Biosample term name'].isin(df_biological_rep)) & subset_intersected['Mapped read length_y'].between(40.0, 200.0)].drop_duplicates(columns)
#        df_2 = combined_experiments
#        combined_experiments = pd.concat([df_2, mapped_experiment_file])

    combined_experiments.to_csv(os.path.join(args.outdir, "metadata.tsv"), sep="\t")

    # save relevant columns into input data lookup for input into ABC code
    prepareLookup(args, combined_experiments, "input_data_lookup")
    prepareLookup(args, metadata_unique, "unique_input_data_lookup")

    return combined_experiments

def downloadFiles(args, df):
    # get download links and start downloading 
    download_links = df[['File download URL_x', 'File download URL_y']].melt(value_name='download_links')
    outfile="linkstodownload.txt"
    download_links[['download_links']].to_csv(os.path.join(args.outdir,  outfile), sep="\t", index=False, header=None)
    
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    commands = "parallel --verbose -j 10 bash downloadFiles.sh ::: {} :::  {}".format(outfile, args.outdir)
    process = subprocess.call(commands, shell=True)
    return outfile 


def get_paired_single_end_datafiles(args, metadata):
    format_files = ["DHS", "H3K27ac", "ATAC"]
    filenames = ["singleend_h3k27ac_dhs_files.tsv", "pairedend_h3k27ac_dhs_files.tsv"]
    for entry in format_files:
        pairedend_data = pd.read_csv(os.path.join(args.outdir, "pairedend_bam_{}_{}.tsv".format(args.genome_assembly, entry)), sep="\t")
        singleend_data = pd.read_csv(os.path.join(args.outdir, "singleend_bam_{}_{}.tsv".format(args.genome_assembly, entry)), sep="\t")
        if entry == "DHS" or "ATAC":
            pairedends_entries = metadata['File accession_x'].loc[metadata['File accession_x'].isin(pairedend_data['File accession'])]
            singleends_entries = metadata['File accession_x'].loc[metadata['File accession_x'].isin(singleend_data['File accession'])]
        if entry == "H3K27ac":
            pairedends_entries = metadata['File accession_y'].loc[metadata['File accession_y'].isin(pairedend_data['File accession'])]
            singleends_entries = metadata['File accession_y'].loc[metadata['File accession_y'].isin(singleend_data['File accession'])]
        
        # write to output files 
        with open(os.path.join(args.outdir, str(filenames[0])), "a") as f:
            for line in singleends_entries:
                f.write(str(line)+".bam")
                f.write("\n")
            f.close()
        with open(os.path.join(args.outdir, str(filenames[1])), "a") as p:
            for line in pairedends_entries:
                p.write(str(line)+".bam")
                p.write("\n")
            p.close()
    grabUnique(args.outdir, filenames)
    

if __name__=="__main__":
    main()
