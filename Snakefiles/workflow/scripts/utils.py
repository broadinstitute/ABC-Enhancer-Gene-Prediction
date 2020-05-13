#! bin/bash 

import sys, os
import argparse 
import pandas as pd
import numpy as np

def is_unique(entries):
    a = entries.to_numpy()
    return (a[0] == a[1:]).all()

def check_x_and_y(entries, column):
    columns = column
    return_values = ['File accession_Accessibility', 'File accession_H3K27ac']
    # check for biological replicates
    if not is_unique(entries[columns[0]]):
        return str(return_values[0]), str(return_values[1])
    elif not is_unique(entries[columns[1]]):
        return str(return_values[1]), str(return_values[0])
    else:
        return -1, -1
    
# This function is used by grabUnique to treat all ambigiously paired bam files as singleend 
def findFilterFiles(outfile, outfile2):
    p = pd.read_csv(outfile, sep="\t", header=None)
    p2 = pd.read_csv(outfile2, sep="\t", header=None)

    common = p[0].loc[p[0].isin(p2[0])]
    remove_files_p2 = p2[0].loc[np.logical_not(p2[0].isin(list(common)))]
    remove_files_p2.to_csv(outfile2, sep="\t", header=False, index=False)

# This function grabs all unique entries in singleend and paireend lookup table for entry into removing duplicates function 
# It also resolves the conflict of bamfiles that are inaccurately labelled as both "pairedend" and "singleend" in ENCODE 
def grabUnique(outdir, filenames):
    for filename in filenames:
        data = pd.read_csv(os.path.join(outdir, filename), sep="\t", header=None)
        unique = data[0].drop_duplicates()
        unique.to_csv(os.path.join(outdir, "unique_"+str(filename)), sep="\t", index=False, header=False)
    # filter both files to ensure that each entry only appears once
    outfile = os.path.join(outdir, "unique_"+str(filenames[0]))
    outfile2 = os.path.join(outdir, "unique_"+str(filenames[1]))

    # find common entries
    findFilterFiles(outfile, outfile2)


# This function prepares the lookup tables for input into ABC
# Where each DHS, H3K27ac bam file is appended to its corresponding celltype
def prepareLookup(args, metadata, title):
    metadata['File accession_Accessibility bam'] = metadata['File accession_Accessibility'].astype('str') + ".nodup.bam"
    metadata['File accession_H3K27ac bam'] = metadata['File accession_H3K27ac'].astype('str') + ".nodup.bam"
    # process biosample term name 
    metadata['Biosample term name'] = [str(i).replace(",", "").replace(" ", "_") for i in metadata['Biosample term name']]
    celltypes = metadata['Biosample term name'].drop_duplicates()
    celltypes.to_csv(os.path.join(args.outdir, "cells.txt"), sep="\t", header=False, index=False)
    metadata[['Biosample term name', 'File accession_Accessibility bam', 'File accession_H3K27ac bam']].to_csv(os.path.join(args.outdir, str(title)+".tsv"), sep="\t", header=False, index=False)
    new_data = metadata.rename(index={key:value for key, value in zip(metadata.index, metadata['Biosample term name'])})
    new_data[['File accession_Accessibility bam', 'File accession_H3K27ac bam']].to_json(os.path.join(args.outdir, "input_data_lookup.json"), orient='index')
    return metadata

def mapExperimentToLength(dhs):
    dhs_lookup = dhs[['Experiment accession', 'Read length', 'Run type']].dropna(subset=['Run type']).drop_duplicates()
    dhs_bam = dhs.loc[(dhs['File format']=='bam') & (dhs['Output type']=='unfiltered alignments')]
    for name, length, paired in zip(dhs_lookup['Experiment accession'], dhs_lookup['Read length'], dhs_lookup['Run type']):
        matched = dhs_bam.loc[dhs_bam['Experiment accession']==name]
        indices = matched.index.astype('int')
        dhs_bam.loc[indices, 'Mapped read length'] = length
        dhs_bam.loc[indices, 'Run type'] = paired
    return dhs_bam
        

# This function updates the lookup table such that entries that have bamfiles with biological replicates undergo a merging process and the collective pooled bam is now updated in the table 
def getExperimentsCombined(metadata, biosample_entries):
    to_combine = {}
    df = metadata
    update_lookup = {}
    for biosample in list(biosample_entries):
        biosample_entry = metadata.loc[metadata['Biosample term name']==biosample]
        if len(biosample_entry) > 1:
            col1, col2 = check_x_and_y(biosample_entry, column=['Biological replicate(s)_Accessibility', 'Biological replicate(s)_H3K27ac'])
            if col1 != -1:
                to_combine[str(biosample).replace(",", "").replace(" ", "_")] = [str(entry)+".nodup.bam" for entry in biosample_entry.loc[:,col1]]
                index = list(biosample_entry.index.astype('int'))
                df.loc[index[0], col1] = str(biosample_entry.loc[index[0], col1])+"_pooled"
                df = df.drop(index[1:])
            else:
                index = list(biosample_entry.index.astype('int'))
                df = df.drop(index[1:])
    return to_combine, df

# This function grabs the samples that have paired ends for paired end bam processing 
# It saves pairedend bams and singleend bams for removal of duplicates 
def obtainDuplicated(args, subset_intersected):
    dhs_duplicates = subset_intersected[subset_intersected.duplicated(['Experiment accession_Accessibility'], keep=False)].drop_duplicates(['Biosample term name', 'Biological replicate(s)_Accessibility'])

    h3k27acduplicates = subset_intersected[subset_intersected.duplicated(['Experiment accession_H3K27ac'], keep=False)].drop_duplicates(['Biosample term name', 'Biological replicate(s)_H3K27ac'])
    
    comb_duplicates = pd.concat([dhs_duplicates, h3k27acduplicates])
    duplicates = comb_duplicates.drop_duplicates()
    duplicates.to_csv(os.path.join(args.outdir, "metadata.tsv"), sep="\t", index=False)
    
    metadata_orig = duplicates.copy()
    # grab celltypes with biological replicates
    df_biological =  duplicates.loc[duplicates.duplicated(['Biosample term name'], keep=False)]
    df_biological.to_csv(os.path.join(args.outdir, "Replicates_metadata.tsv"), sep="\t", index=False)
    df_biological_rep = df_biological['Biosample term name'].drop_duplicates()
    to_combine, metadata_unique = getExperimentsCombined(duplicates, df_biological_rep)
    with open(os.path.join(args.outdir, "Experiments_ToCombine.txt.tmp"), "w") as f:
        for key, value in zip(to_combine.keys(), to_combine.values()):
            f.write(str(key))
            f.write("\t")
            f.write(str(list(value)))
            f.write("\n")
        f.close()
    
    return metadata_orig, metadata_unique

