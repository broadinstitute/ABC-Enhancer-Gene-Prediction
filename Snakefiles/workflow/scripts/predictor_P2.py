import numpy as np
import pandas as pd
from tools import *
import sys, os
import time
import pyranges as pr
from hic import *

def make_predictions(chromosome, enhancers, genes, args):
    pred = make_pred_table(chromosome, enhancers, genes, args)
    pred = annotate_predictions(pred, args.tss_slop)
    pred = add_powerlaw_to_predictions(pred, args)

    #if Hi-C directory is not provided, only powerlaw model will be computed
    if args.HiCdir:
        hic_file, hic_norm_file, hic_is_vc = get_hic_file(chromosome, args.HiCdir, hic_type = args.hic_type)
        pred = add_hic_to_enh_gene_table(enhancers, genes, pred, hic_file, hic_norm_file, hic_is_vc, chromosome, args)
        pred = compute_score(pred, [pred['activity_base'], pred['hic_contact_pl_scaled_adj']], "ABC")
        pred = compute_synergy_score(pred, [pred['activity_base'], pred['hic_contact_pl_scaled_adj']], "ABC")

    pred = compute_score(pred, [pred['activity_base'], pred['powerlaw_contact_reference']], "powerlaw")
    #pred = compute_synergy_score(pred, [pred['activity_base'], pred['powerlaw_contact_reference']], "powerlaw")
    print("N row after making predictions: " + str(len(pred)))
    return pred

def make_pred_table(chromosome, enh, genes, args):
    print('Making putative predictions table...')
    t = time.time()
    enh['enh_midpoint'] = (enh['start'] + enh['end'])/2
    enh['enh_idx'] = enh.index
    genes['gene_idx'] = genes.index
    enh_pr = df_to_pyranges(enh)

    genes_pr = df_to_pyranges(genes, start_col = 'TargetGeneTSS', end_col = 'TargetGeneTSS', start_slop=args.window, end_slop = args.window)

    pred = enh_pr.join(genes_pr).df.drop(['Start_b','End_b','chr_b','Chromosome','Start','End'], axis = 1)
    pred['distance'] = abs(pred['enh_midpoint'] - pred['TargetGeneTSS'])
    pred = pred.loc[pred['distance'] < args.window,:] #for backwards compatability
    #without pyranges version
    # else:
    #     enh['temp_merge_key'] = 0
    #     genes['temp_merge_key'] = 0

    #     #Make cartesian product and then subset to EG pairs within window. 
    #     #TO DO: Replace with pyranges equivalent of bedtools intersect or GRanges overlaps 
    #     pred = pd.merge(enh, genes, on = 'temp_merge_key')

    #     pred['enh_midpoint'] = (pred['start'] + pred['end'])/2
    #     pred['distance'] = abs(pred['enh_midpoint'] - pred['TargetGeneTSS'])
    #     pred = pred.loc[pred['distance'] < args.window,:]

    #     print('Done. There are {} putative enhancers for chromosome {}'.format(pred.shape[0], chromosome))
    #     print('Elapsed time: {}'.format(time.time() - t))

    return pred

def add_hic_to_enh_gene_table(enh, genes, pred, hic_file, hic_norm_file, hic_is_vc, chromosome, args):
    print('Begin HiC')
    HiC = load_hic(hic_file = hic_file, 
                    hic_norm_file = hic_norm_file,
                    hic_is_vc = hic_is_vc,
                    hic_type = args.hic_type, 
                    hic_resolution = args.hic_resolution, 
                    tss_hic_contribution = args.tss_hic_contribution, 
                    window = args.window, 
                    min_window = 0, 
                    gamma = args.hic_gamma,
                    scale = args.hic_scale)

    #Add hic to pred table
    #At this point we have a table where each row is an enhancer/gene pair. 
    #We need to add the corresponding HiC matrix entry.
    #If the HiC is provided in juicebox format (ie constant resolution), then we can just merge using the indices
    #But more generally we do not want to assume constant resolution. In this case hic should be provided in bedpe format

    t = time.time()
    if args.hic_type == "bedpe":
        #Use pyranges to compute overlaps between enhancers/genes and hic bedpe table
        #Consider each range of the hic matrix separately - and merge each range into both enhancers and genes. 
        #Then remerge on hic index

        HiC['hic_idx'] = HiC.index
        hic1 = df_to_pyranges(HiC, start_col='x1', end_col='x2', chr_col='chr1')
        hic2 = df_to_pyranges(HiC, start_col='y1', end_col='y2', chr_col='chr2')

        #Overlap in one direction
        enh_hic1 = df_to_pyranges(enh, start_col = 'enh_midpoint', end_col = 'enh_midpoint', end_slop = 1).join(hic1).df
        genes_hic2 = df_to_pyranges(genes, start_col = 'TargetGeneTSS', end_col = 'TargetGeneTSS', end_slop = 1).join(hic2).df
        ovl12 = enh_hic1[['enh_idx','hic_idx','hic_contact']].merge(genes_hic2[['gene_idx', 'hic_idx']], on = 'hic_idx')

        #Overlap in the other direction
        enh_hic2 = df_to_pyranges(enh, start_col = 'enh_midpoint', end_col = 'enh_midpoint', end_slop = 1).join(hic2).df
        genes_hic1 = df_to_pyranges(genes, start_col = 'TargetGeneTSS', end_col = 'TargetGeneTSS', end_slop = 1).join(hic1).df
        ovl21 = enh_hic2[['enh_idx','hic_idx','hic_contact']].merge(genes_hic1[['gene_idx', 'hic_idx']], on = ['hic_idx'])

        #Concatenate both directions and merge into preditions
        ovl = pd.concat([ovl12, ovl21]).drop_duplicates()
        pred = pred.merge(ovl, on = ['enh_idx', 'gene_idx'], how = 'left')
        pred.fillna(value={'hic_contact' : 0}, inplace=True)
    elif args.hic_type == "juicebox" or args.hic_type == "avg" or args.hic_type == "avg_hg19":
        #Merge directly using indices
        #Could also do this by indexing into the sparse matrix (instead of merge) but this seems to be slower
        #Index into sparse matrix
        #pred['hic_contact'] = [HiC[i,j] for (i,j) in pred[['enh_bin','tss_bin']].values.tolist()]
        
        pred['enh_bin'] = np.floor(pred['enh_midpoint'] / args.hic_resolution).astype(int)
        pred['tss_bin'] = np.floor(pred['TargetGeneTSS'] / args.hic_resolution).astype(int)
        if not hic_is_vc:
            #in this case the matrix is upper triangular.
            #
            pred['bin1'] = np.amin(pred[['enh_bin', 'tss_bin']], axis = 1)
            pred['bin2'] = np.amax(pred[['enh_bin', 'tss_bin']], axis = 1)
            pred = pred.merge(HiC, how = 'left', on = ['bin1','bin2'])
            pred.fillna(value={'hic_contact' : 0}, inplace=True)
        else:
            # The matrix is not triangular, its full
            # For VC assume genes correspond to rows and columns to enhancers
            pred = pred.merge(HiC, how = 'left', left_on = ['tss_bin','enh_bin'], right_on=['bin1','bin2'])

        pred.fillna(value={'hic_contact' : 0}, inplace=True)

        # QC juicebox HiC
        pred = qc_hic(pred)

    pred.drop(['x1','x2','y1','y2','bin1','bin2','enh_idx','gene_idx','hic_idx','enh_midpoint','tss_bin','enh_bin'], inplace=True, axis = 1, errors='ignore')
        
    print('HiC added to predictions table. Elapsed time: {}'.format(time.time() - t))
    # Add powerlaw scaling
    pred = scale_hic_with_powerlaw(pred, args)

    #Add pseudocount
    pred = add_hic_pseudocount(pred, args)

    print("HiC Complete")
    #print('Elapsed time: {}'.format(time.time() - t))

    return(pred)

def scale_hic_with_powerlaw(pred, args):
    #Scale hic values to reference powerlaw

    if not args.scale_hic_using_powerlaw:
#        values = pred.loc[pred['hic_contact']==0].index.astype('int')
#        pred.loc[values, 'hic_contact'] = pred.loc[values, 'powerlaw_contact']
        pred['hic_contact_pl_scaled'] = pred['hic_contact']
    else:
        pred['hic_contact_pl_scaled'] = pred['hic_contact'] * (pred['powerlaw_contact_reference'] / pred['powerlaw_contact'])

    return(pred)

def add_powerlaw_to_predictions(pred, args):
    pred['powerlaw_contact'] = get_powerlaw_at_distance(pred['distance'].values, args.hic_gamma, args.hic_scale)
    pred['powerlaw_contact_reference'] = get_powerlaw_at_distance(pred['distance'].values, args.hic_gamma_reference, args.hic_scale)

    return pred

def add_hic_pseudocount(pred, args):
    # Add a pseudocount based on the powerlaw expected count at a given distance

    powerlaw_fit = get_powerlaw_at_distance(pred['distance'].values, args.hic_gamma)
    powerlaw_fit_at_ref = get_powerlaw_at_distance(args.hic_pseudocount_distance, args.hic_gamma)
    
    pseudocount = np.amin(pd.DataFrame({'a' : powerlaw_fit, 'b' : powerlaw_fit_at_ref}), axis = 1)
    pred['hic_pseudocount'] = pseudocount
    pred['hic_contact_pl_scaled_adj'] = pred['hic_contact_pl_scaled'] + pseudocount

    return(pred)

def qc_hic(pred, threshold = .01):
    # Genes with insufficient hic coverage should get nan'd

    summ = pred.loc[pred['isSelfPromoter'],:].groupby(['TargetGene']).agg({'hic_contact' : 'sum'})
    bad_genes = summ.loc[summ['hic_contact'] < threshold,:].index

    pred.loc[pred['TargetGene'].isin(bad_genes), 'hic_contact'] = np.nan

    return pred

def compute_score(enhancers, product_terms, prefix):

    scores = np.column_stack(product_terms).prod(axis = 1)

    enhancers[prefix + '.Score.Numerator'] = scores
    enhancers[prefix + '.Score'] = enhancers[prefix + '.Score.Numerator'] / enhancers.groupby(['TargetGene','TargetGeneTSS'])[prefix + '.Score.Numerator'].transform('sum')
    
    pd.set_option('display.max_columns', None)
    #print(enhancers.head(10))
    return(enhancers)


def compute_synergy_score(enhancers, product_terms, prefix):
    window = 5e6
    P2_coefficient = 0.708876033750941
    if prefix=='ABC':
        contact_col = 'hic_contact_pl_scaled_adj'
    else:
        contact_col = 'powerlaw_contact_reference'

    # calculate wt gex (same as denominator of ABC score)
    enhancers[prefix + '.wt.gex'] = enhancers.groupby(['TargetGene','TargetGeneTSS'])[prefix + '.Score.Numerator'].transform('sum')
    enhancers['enh_midpoint'] = (enhancers['start'] + enhancers['end'])/2

    # make E-E matrix
    t = time.time()
    print(enhancers.columns)
    print('Computing enhancer-enhancer matrix...')
    EE_matrix = compute_EE_matrix(enhancers)
    print('Computed enhancer-enhancer matrix. Elapsed time: {}'.format(time.time() - t))

    # unique gene/TSS pairs
    genes = enhancers[['TargetGene', 'TargetGeneTSS']].drop_duplicates().reset_index(drop=True)
    print('N genes: ' + str(len(genes)))
    #print(enhancers.head())
    all_pred_list = []

    # loop through each gene
    for i in range(len(genes)):
        this_gene = genes.loc[i, 'TargetGene']
        #print('this_gene: ' + this_gene)
        this_TSS = genes.loc[i, 'TargetGeneTSS']
        #print('this_TSS: ' + str(this_TSS))

    # testing
    #for i in range(1):
        #this_gene = 'SAMD11'
        #this_TSS = 925740

        # subset predictions to this gene
        this_pred = enhancers.loc[(enhancers['TargetGene']==this_gene) & (enhancers['TargetGeneTSS']==this_TSS)].reset_index(drop=True)
        #print('N row this_pred: ' + str(len(this_pred)))
        # subset enhancer matrix to enhancers where Target and Element are within 5 Mb of TSS
        this_EE = EE_matrix.loc[(np.abs(EE_matrix['enh_midpoint_target']-this_TSS)<window) & (np.abs(EE_matrix['enh_midpoint_element']-this_TSS)<window)] 
        #print('N row this_EE: ' + str(len(this_EE)))
        # add element-gene contact to EE_matrix
        this_pred['start_element'] = this_pred['start']
        this_pred['end_element'] = this_pred['end']
        # set noselfprom contact to 0 if element is a self promoter
        this_pred[prefix + '.contact.noselfprom'] = [0 if x==True else y for x,y in zip(this_pred['isSelfPromoter'],this_pred[contact_col])]

        EE_contact = pd.merge(this_EE, this_pred[['start_element', 'end_element', contact_col, prefix+'.contact.noselfprom']], how="left", on=['start_element', 'end_element'])
        #print('N row EE_contact: ' + str(len(EE_contact)))
        # loop thorugh each enhancer in this_pred 
        this_pred[prefix + '.ko.gex'] = 0
        this_pred[prefix + '.ko.gex.noselfprom'] = 0
        this_pred[prefix + '.wt.gex.noselfprom'] = 0

        for j in range(len(this_pred)):
            this_start = this_pred.loc[j, 'start']
            this_end = this_pred.loc[j, 'end']
            EE_subset_file = '/oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/2023_0201_Predictions_test/EE_subset_' + this_gene + '_' + str(this_TSS) + '_' + str(this_start) + '_' + str(this_end) + '.tsv'
            EE_subset = EE_contact.loc[(EE_contact['start_target']==this_start) & (EE_contact['end_target']==this_end)]
            this_pred.loc[j, prefix + '.ko.gex'] = np.sum(np.multiply(EE_subset['ko_activity'], EE_subset[contact_col]))
            this_pred.loc[j, prefix + '.ko.gex.noselfprom'] = np.sum(np.multiply(EE_subset['ko_activity'], EE_subset[prefix + '.contact.noselfprom']))
            this_pred.loc[j, prefix + '.wt.gex.noselfprom'] = np.sum(np.multiply(EE_subset['activity_base_element'], EE_subset[prefix + '.contact.noselfprom']))

            #EE_subset.to_csv(EE_subset_file, sep="\t", index=False, header=True, float_format="%.6f", na_rep="NaN")

        this_pred = this_pred.drop(columns=['start_element', 'end_element'])
        all_pred_list.append(this_pred)
    all_pred = pd.concat(all_pred_list)
    # calculate synergy score
    all_pred[prefix + '.synergy.Score'] = (all_pred[prefix + '.wt.gex'] - all_pred[prefix + '.ko.gex'])/all_pred[prefix + '.wt.gex']
    all_pred[prefix + '.synergy.episelfprom.Score'] = (all_pred[prefix + '.wt.gex.noselfprom'] - all_pred[prefix + '.ko.gex.noselfprom'])/all_pred[prefix + '.wt.gex.noselfprom']

    all_pred[prefix + '.P2.synergy.Score'] = P2_coefficient ** all_pred['P2_indicator'] * all_pred[prefix + '.synergy.Score']

    all_pred[prefix + '.synergy-normal.diff'] = all_pred[prefix + '.synergy.Score'] - all_pred[prefix + '.Score']

    return(all_pred)

# compute enhancer-enhancer contact and KO activity for all enhancers within 5Mb of target enhancer
def compute_EE_matrix(enhancers):
    window = 5e6
    gamma = 1.0044322152129266
    x0 = 5000
    scale = 5.522049023040865
    b = 0.733487873281097
    Ci_efficiency = 0.8

    pl_contact = lambda dist: np.exp(scale + -1*gamma * np.log(max(x0, dist))) / np.exp(scale + -1*gamma * np.log(x0))

    # get list of unique enhancers
    print('N predictions: ' + str(len(enhancers)))
    enhancers = enhancers[['start', 'end', 'activity_base', 'enh_midpoint', 'isPromoterNotP2']].drop_duplicates().reset_index(drop=True)
    print('N enhancers: ' + str(len(enhancers)))

    all_enh_list = []
    # loop through enhancers 
    for i in range(len(enhancers)):
        this_enh = enhancers.iloc[i]
        # pull all enhancers within [window] of target enhancer
        this_enh_matrix = enhancers.loc[(abs(this_enh['enh_midpoint']-enhancers['enh_midpoint'])<=window)]
        this_enh_matrix = this_enh_matrix[['start', 'end', 'activity_base']]
        # name columns for location & activity of each enhancer in pair
        isPromoterNotP2_target = this_enh['isPromoterNotP2']
        this_enh_matrix = this_enh_matrix.assign(start_target=this_enh['start'])
        this_enh_matrix['start_element'] = this_enh_matrix['start']
        this_enh_matrix.loc[:, 'end_target'] = this_enh['end']
        this_enh_matrix['end_element'] = this_enh_matrix['end']
        this_enh_matrix['enh_midpoint_target'] = (this_enh_matrix['start_target'] + this_enh_matrix['end_target'])/2
        this_enh_matrix['enh_midpoint_element'] = (this_enh_matrix['start_element'] + this_enh_matrix['end_element'])/2
        this_enh_matrix['activity_base_target'] = this_enh['activity_base']
        this_enh_matrix['activity_base_element'] = this_enh_matrix['activity_base']
        # calculate distance and contact between enhancers
        this_enh_matrix['enh-enh_distance'] = np.abs(this_enh_matrix['enh_midpoint_target']-this_enh_matrix['enh_midpoint_element'])
        this_enh_matrix['enh-enh_contact'] = [pl_contact(dist) for dist in this_enh_matrix['enh-enh_distance']] 
        # calculate KO activity at each element given target is perturbed (min = 0)
        if isPromoterNotP2_target: # if target is a promoter for a not P2 gene, perturbing it doesn't affect anything
            this_enh_matrix['ko_activity'] = this_enh_matrix['activity_base_element']
        else:
            this_enh_matrix['ko_activity'] = this_enh_matrix['activity_base_element'] - b * Ci_efficiency * np.multiply(this_enh_matrix['enh-enh_contact'], this_enh_matrix['activity_base_target'])
        this_enh_matrix['ko_activity'] = this_enh_matrix['ko_activity'].clip(0, None)

        all_enh_list.append(this_enh_matrix)
    
    all_enh = pd.concat(all_enh_list)
    return(all_enh)


def annotate_predictions(pred, tss_slop=500):
    #TO DO: Add is self genic
    pred['isSelfPromoter'] = np.logical_and.reduce((pred['class'] == 'promoter' , pred.start - tss_slop < pred.TargetGeneTSS, pred.end + tss_slop > pred.TargetGeneTSS))

    return(pred)

def make_gene_prediction_stats(pred, args):
    summ1 = pred.groupby(['chr','TargetGene','TargetGeneTSS']).agg({'TargetGeneIsExpressed' : lambda x: set(x).pop(), args.score_column : lambda x: all(np.isnan(x)) ,  'name' : 'count'})
    summ1.columns = ['geneIsExpressed', 'geneFailed','nEnhancersConsidered']

    summ2 = pred.loc[pred['class'] != 'promoter',:].groupby(['chr','TargetGene','TargetGeneTSS']).agg({args.score_column : lambda x: sum(x > args.threshold)})
    summ2.columns = ['nDistalEnhancersPredicted']
    summ1 = summ1.merge(summ2, left_index=True, right_index=True)

    summ1.to_csv(os.path.join(args.outdir, "GenePredictionStats.txt"), sep="\t", index=True)
