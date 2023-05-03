import pandas as pd
import numpy as np
import os, sys

## Specify input files here 
# K562_ID_2644
feature=sys.argv[1]
avg_track="02212022_hg38Preds_QNORM-avg_hic_track1"
refgenes = "/oak/stanford/groups/akundaje/kmualim/ABC-Enhancer-Gene-Prediction/reference/hg38/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed"
predictionfile = "/oak/stanford/groups/akundaje/kmualim/02212022_hg38Preds_QNORM-hic/{}/Predictions/EnhancerPredictionsAllPutative.txt.gz".format(feature)
enhoutput = "/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/EGFeatureList.txt".format(avg_track,feature)
filename="/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/EnhancerRegions_extended_RefSeq.intersected.txt.gz".format(avg_track,feature)
extended_enhregions = "/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/EnhancerRegions_extended.txt".format(avg_track,feature)
enhancerlist = "/oak/stanford/groups/akundaje/kmualim/{}/{}/Neighborhoods/EnhancerList.txt".format(avg_track,feature)
candidate_enhancersfile = "/oak/stanford/groups/akundaje/kmualim/{}/{}/Neighborhoods/EnhancerList_noPromoter.txt".format(avg_track,feature)
filename_enh = "/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/EnhancerRegions_extended_CandidateReg.txt.gz".format(avg_track,feature)

### Usage: python generate_feature_list.py /oak/stanford/groups/akundaje/kmualim/02212022_hg38Preds_QNORM-hic/K562_ID_2644/Predictions/EnhancerPredictionsAllPutative.txt.gz /oak/stanford/groups/akundaje/kmualim/02212022_hg38Preds_QNORM-hic/K562_ID_2644/Predictions/EnhancerRegions_extended.txt /oak/stanford/groups/akundaje/kmualim/02212022_hg38Preds_QNORM-hic/K562_ID_2644/Predictions/K562_EGFeatureList.txt ../../references/hg38/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed /oak/stanford/groups/akundaje/kmualim/02212022_hg38Preds_QNORM-hic/K562_ID_2644/Predictions/EnhancerRegions_extended_RefSeq.intersected.txt  

# Read in input prediction file
subset = pd.read_csv(enhancerlist, sep="\t", usecols=['chr', 'start', 'end', 'name'])
print(len(subset))
##### columns in input prediction file: chr	start	end	name	class	activity_base	TargetGene	TargetGeneTSS	TargetGeneExpression	TargetGenePromoterActivityQuantile	TargetGeneIsExpressed	distance	isSelfPromoter	powerlaw_contact	powerlaw_contact_reference	hic_contact	hic_contact_pl_scaled	hic_pseudocount	hic_contact_pl_scaled_adj	ABC.Score.Numerator	ABC.Score	powerlaw.Score.Numerator	powerlaw.Score	CellType
#### subsetting data to just columns required
#subset = data[['chr', 'start', 'end', 'name', 'class', 'activity_base', 'TargetGene', 'TargetGeneTSS', 'TargetGeneExpression', 'TargetGenePromoterActivityQuantile', 'TargetGeneIsExpressed', 'distance', 'hic_contact', 'powerlaw_contact', 'hic_contact_pl_scaled', 'ABC.Score.Numerator', 'ABC.Score', 'CellType']]
##### remove entries that are promoters
subset = subset.loc[subset['class']!="promoter"]
if len(subset) > 1:
    ##### Get midpoint of enhancer region
    subset['midpoint'] = subset['start']+0.5*(subset['end']-subset['start'])
    subset['midpoint'] = subset['midpoint'].astype('int')
    
    ##### Make the end be midpoint of enhancer + distance (This gives you the end coordinate of distance range)
    subset['end1'] = subset['midpoint']+subset['distance']
    
    ## If gene is located upstream of enhancer, modify the start to be the beginning of the TargetGeneTSS and the end be the midpoint of the enhancer
    entries = subset.loc[subset['TargetGeneTSS'] < subset['midpoint']].index.astype('int')
    subset.loc[entries, 'end1'] = subset.loc[entries, 'end']
    subset.loc[entries, 'start'] = subset.loc[entries, 'TargetGeneTSS']
    subset['start'] = subset['start'].astype('int')
    subset['end1'] = subset['end1'].astype('int')
    # This file will be used to intersect gene annotations to count How many protein-coding TSSs away is the enhancer from the promoter?  (i.e., how many protein-coding gene TSSs are located between the enhancer and promoter?  0 = closest TSS)
    # Intersect this file with gene TSS
    subset[['chr', 'start', 'end1', 'name', 'TargetGene']].to_csv(extended_enhregions, sep="\t", index=False)
    
    ## Intersect midpoint of enhancer to target gene regions with GeneTSS
    ## This will include overlaps with the TargetGene 
    os.system("sed '1d' {} | bedtools intersect -a stdin -b {} -wa -wb | cut -f4,5 | gzip > {}".format(extended_enhregions, refgenes,  filename))
    ##
    candidate_enhancers = pd.read_csv(enhancerlist, sep="\t")
    candidate_enhancers = candidate_enhancers.loc[candidate_enhancers['class']!="promoter"]
    candidate_enhancers[['chr', 'start', 'end']].to_csv(candidate_enhancersfile, sep="\t", index=False)
    #### How many candidate enhancer regions are located between the enhancer and the promoter? 
    os.system("sed '1d' {} | bedtools intersect -a stdin -b {} -wa -wb | cut -f4,5 | gzip > {}".format(extended_enhregions, candidate_enhancersfile,  filename_enh))
    print("Reading in {}".format(filename_enh))
    predictions_can = pd.read_csv(filename_enh, sep="\t", names=['class', 'gene'])
    num_candidate_enhancers = predictions_can.groupby(['class', 'gene']).size().reset_index()
    num_candidate_enhancers.to_csv("/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/NumCandidateEnhGene.txt".format(avg_track, feature), sep="\t", index=False)
    print("Saved num candidate enhancers")
    
    print("Reading in {}".format(filename))
    ## Read in intersected midpoint of enhancer to target gene regions with GeneTSS file
    predictions = pd.read_csv(filename, sep="\t", names=['class', 'gene'])
    # Calculate the number of TSS regions that fall within the enhancer to target gene regions. 
    num_tss_between_enh_and_gene = predictions.groupby(['class', 'gene']).size().reset_index()
    num_tss_between_enh_and_gene.to_csv("/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/NumTSSEnhGene.txt".format(avg_track, feature), sep="\t", index=False)
    print("Saved num TSS between enh and gene")

    ############ Generate Num/Sum Enhancers within 5kb/10kb ############
    ##### Make the end be midpoint of enhancer + distance (This gives you the end coordinate of distance range)
    subset[['chr', 'midpoint', 'midpoint', 'name']].to_csv("/oak/stanford/groups/akundaje/kmualim/{}/{}/Neighborhoods/EnhancerList_midpoint.bed".format(avg_track, feature), sep="\t", index=False, header=False)
    os.system("bedtools slop -b 5000 -i /oak/stanford/groups/akundaje/kmualim/{}/{}/Neighborhoods/EnhancerList_midpoint.bed -g /oak/stanford/groups/akundaje/kmualim/ABC-Enhancer-Gene-Prediction/sherlock_scripts/hg38.chrom.sizes > /oak/stanford/groups/akundaje/kmualim/{}/{}/Neighborhoods/EnhancerRegions.slop_5kb.bed".format(avg_track, feature, avg_track, feature))
    os.system("bedtools slop -b 10000 -i /oak/stanford/groups/akundaje/kmualim/{}/{}/Neighborhoods/EnhancerList_midpoint.bed -g /oak/stanford/groups/akundaje/kmualim/ABC-Enhancer-Gene-Prediction/sherlock_scripts/hg38.chrom.sizes > /oak/stanford/groups/akundaje/kmualim/{}/{}/Neighborhoods/EnhancerRegions.slop_10kb.bed".format(avg_track, feature, avg_track, feature))
    os.system("zcat {} | cut -f1,2,3,4,6 | sed '1d' > /oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/EnhancerPredictionsAllPutative.slim.bed".format(predictionfile, avg_track, feature))
    os.system("bedtools intersect -a /oak/stanford/groups/akundaje/kmualim/{}/{}/Neighborhoods/EnhancerRegions.slop_5kb.bed -b /oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/EnhancerPredictionsAllPutative.slim.bed -wa -wb > /oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/NumEnhancers5kb.txt".format(avg_track, feature, avg_track, feature, avg_track, feature))
    os.system("bedtools intersect -a /oak/stanford/groups/akundaje/kmualim/{}/{}/Neighborhoods/EnhancerRegions.slop_10kb.bed -b /oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/EnhancerPredictionsAllPutative.slim.bed -wa -wb > /oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/NumEnhancers10kb.txt".format(avg_track, feature, avg_track, feature, avg_track, feature))
    data = pd.read_csv("/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/NumEnhancers5kb.txt".format(avg_track, feature), sep="\t", header=None)
    data1 = data.loc[data[3]!=data[7]]
    data2 = data1.groupby([3]).size().reset_index(name='count')
    data3 = data1.groupby([3])[8].sum().reset_index(name='sum')
    data2.to_csv("/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/NumEnhancersEG5kb.txt".format(avg_track, feature), sep="\t", header=False, index=False)
    data3.to_csv("/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/SumEnhancersEG5kb.txt".format(avg_track, feature), sep="\t", header=False, index=False)
    data = pd.read_csv("/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/NumEnhancers10kb.txt".format(avg_track, feature), sep="\t", header=None)
    data1 = data.loc[data[3]!=data[7]]
    data2 = data1.groupby([3]).size().reset_index(name='count')
    data3 = data1.groupby([3])[8].sum().reset_index(name='sum')
    data2.to_csv("/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/NumEnhancersEG10kb.txt".format(avg_track, feature), sep="\t", header=False, index=False)
    data3.to_csv("/oak/stanford/groups/akundaje/kmualim/{}/{}/Predictions/SumEnhancersEG10kb.txt".format(avg_track, feature), sep="\t", header=False, index=False)

