import os,sys
import pandas as pd
import seaborn as sns
import glob
import numpy as np
from subprocess import check_call, check_output, PIPE, Popen, getoutput, CalledProcessError

def grabStatistics(arr_values):
    mean = arr_values.mean()
    median = arr_values.median()
    std = arr_values.std()
    return mean, median, std

# Generates QC Prediction Metrics:
def GrabQCMetrics(prediction_df, outdir):
    EnhancerPerGene = prediction_df.groupby(['TargetGene']).size()
    EnhancerPerGene.to_csv(os.path.join(outdir,"EnhancerPerGene.txt"), sep="\t")

    # Grab Number of Enhancers Per Gene
    GeneMean, GeneMedian, GeneStdev = grabStatistics(EnhancerPerGene)

    # Grab Number of genes per enhancers
    NumGenesPerEnhancer = prediction_df[['chr', 'start', 'end']].groupby(['chr', 'start', 'end']).size()
    NumGenesPerEnhancer.to_csv(os.path.join(outdir,"GenesPerEnhancer.txt"), sep="\t")
    mean_genes_per_enhancer, median_genes_per_enhancer,stdev_genes_per_enhancer = grabStatistics(NumGenesPerEnhancer) 

    # Grab Number of Enhancer-Gene Pairs Per Chromsome
    enhancergeneperchrom = prediction_df.groupby(['chr']).size()
    enhancergeneperchrom.to_csv(os.path.join(outdir, "EnhancerGenePairsPerChrom.txt"), sep="\t")
    
    mean_enhancergeneperchrom, median_enhancergeneperchrom ,stdev_enhancergeneperchrom = grabStatistics(enhancergeneperchrom) 

    # Enhancer-Gene Distancee
    distance = np.array(prediction_df['distance'])
    thquantile = np.percentile(distance, 10)
    testthquantile = np.percentile(distance, 90)

    
    # Plot Distributions and save as png
    PlotDistribution(NumGenesPerEnhancer, "NumberOfGenesPerEnhancer", outdir)
    PlotDistribution(EnhancerPerGene, "NumberOfEnhancersPerGene", outdir)
    PlotDistribution(enhancergeneperchrom, "EnhancersPerChromosome", outdir)
    PlotDistribution(distance, "EnhancerGeneDistance", outdir)

    pred_metrics={}
    pred_metrics['MedianEnhPerGene'] = GeneMedian
    pred_metrics['StdEnhPerGene'] = GeneStdev
    pred_metrics['MedianGenePerEnh'] = median_genes_per_enhancer
    pred_metrics['StdGenePerEnh'] = stdev_genes_per_enhancer
    pred_metrics['MeanEnhPerGene'] = GeneMean
    pred_metrics['MeanGenePerEnh'] = mean_genes_per_enhancer
    pred_metrics['MedianEGDist'] = np.median(distance)
    pred_metrics['MeanEGDist'] = np.mean(distance)
    pred_metrics['StdEGDist'] = np.std(distance)
    pred_metrics['NumEnhancersPerChrom'] = median_enhancergeneperchrom
    pred_metrics['MeanEnhancersPerChrom'] = mean_enhancergeneperchrom
    pred_metrics['EG10th'] = thquantile
    pred_metrics['EG90th'] = testthquantile
    return pred_metrics

def PlotQuantilePlot(EnhancerList, title, outdir):
    i='DHS'
    ax = sns.scatterplot('DHS.RPM', 'DHS.RPM.quantile', data=EnhancerList)
    ax.set_title(title)
    ax.set_ylabel('RPM.quantile')
    ax.set_xlabel('RPM')
    fig = ax.get_figure()
    outfile = os.path.join(outdir, i+str(title)+".pdf")
    fig.savefig(outfile, format='pdf')
    
    i="H3K27ac"
    ax = sns.scatterplot('H3K27ac.RPM', 'H3K27ac.RPM.quantile', data=EnhancerList)
    ax.set_title(title)
    ax.set_ylabel('RPM.quantile')
    ax.set_xlabel('RPM')
    fig = ax.get_figure()
    outfile = os.path.join(outdir, i+str(title)+".pdf")
    fig.savefig(outfile, format='pdf')

def NeighborhoodFileQC(pred_metrics, neighborhood_dir, outdir, feature):
    x = glob.glob(os.path.join(neighborhood_dir, "Enhancers.{}.*CountReads.bedgraph".format(feature)))
    y = glob.glob(os.path.join(neighborhood_dir, "Genes.TSS1kb.{}.*CountReads.bedgraph".format(feature)))
    z = glob.glob(os.path.join(neighborhood_dir, "Genes.{}.*CountReads.bedgraph".format(feature)))
    
    data = pd.read_csv(x[0],sep="\t", header=None)
    data1 = pd.read_csv(y[0],sep="\t", header=None)
    data2 = pd.read_csv(z[0],sep="\t", header=None)

    counts = data.iloc[:, 3].sum()
    counts2 = data1.iloc[:,3].sum()
    counts3 = data2.iloc[:,3].sum()
   
    pred_metrics['countsEnhancers_{}'.format(feature)] = counts
    pred_metrics['countsGeneTSS_{}'.format(feature)] = counts2
    pred_metrics['countsGenes_{}'.format(feature)] = counts3
    return pred_metrics

# Generates peak file metrics
def PeakFileQC(pred_metrics, macs_peaks, outdir):
    if macs_peaks.endswith(".gz"):
        peaks = pd.read_csv(macs_peaks, compression="gzip", sep="\t", header=None)
    else:
        peaks = pd.read_csv(macs_peaks, sep="\t", header=None)

    # Calculate metrics for candidate regions 
    outfile = os.path.join(outdir, os.path.basename(macs_peaks) + ".sorted.candidateRegions.bed")
    
    candidateRegions = pd.read_csv(outfile, sep="\t", header=None)
    candidateRegions['dist'] = candidateRegions[2] - candidateRegions[1]
    candreg = list(candidateRegions['dist'])
    PlotDistribution(candreg, 'WidthOfCandidateRegions', outdir)
    

    # Calculate width of peaks 
    peaks['dist'] = peaks[2]-peaks[1]
    peaks_array = list(peaks['dist'])
    pred_metrics['NumPeaks'] = len(peaks['dist'])
    pred_metrics['MedWidth'] = peaks['dist'].median()
    pred_metrics['MeanWidth'] = peaks['dist'].mean()
    pred_metrics['StdWidth'] = peaks['dist'].std()
    pred_metrics['NumCandidate'] = len(candidateRegions['dist'])
    pred_metrics['MedWidthCandidate'] = candidateRegions['dist'].median()
    pred_metrics['MeanWidthCandidate'] = candidateRegions['dist'].mean()
    pred_metrics['StdWidthCandidate'] = candidateRegions['dist'].std()

    return pred_metrics

# Plots and saves a distribution as *.png
def PlotDistribution(array, title, outdir):
    ax = sns.distplot(array)
    ax.set_title(title)
    ax.set_ylabel('Estimated PDF of distribution')
    ax.set_xlabel('Counts')
    fig = ax.get_figure()
    outfile = os.path.join(outdir, str(title)+".pdf")
    fig.savefig(outfile, format='pdf')
