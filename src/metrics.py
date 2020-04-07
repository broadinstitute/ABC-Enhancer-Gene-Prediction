import os,sys
import pandas as pd
import seaborn as sns
import glob
import numpy as np
from subprocess import check_call, check_output, PIPE, Popen, getoutput, CalledProcessError

def grab_nearest_tss_from_peak(macs_peaks, genome_tss, outdir):
    # Grab nearest tss from peak
    outfile = os.path.join(outdir, os.path.basename(macs_peaks))
    files = pd.read_csv(outfile, sep="\t")
    annotated_peaks = os.path.join(outdir, os.path.basename(macs_peaks) + ".annotated_peaks.bed")
    if outfile.endswith(".gz"):
        sort_command = "zcat {outfile} | sort -k1,1 -k2,2n > {outfile}.sorted"
    else:
        sort_command = "sort -k1,1 -k2,2n {outfile} > {outfile}.sorted"
    sort_command = sort_command.format(**locals())
    p = Popen(sort_command, stdout=PIPE, stderr=PIPE, shell=True)
    print("Running:" + sort_command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    sort_genome_command = "sort -k1,1 -k2,2n {genome_tss} > {genome_tss}.sorted"
    sort_genome_command = sort_genome_command.format(**locals())
    p = Popen(sort_genome_command, stdout=PIPE, stderr=PIPE, shell=True)
    print("Running:" + sort_genome_command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')

    command = "bedtools closest -a {outfile}.sorted -b {genome_tss}.sorted -d > {annotated_peaks}"
    command = command.format(**locals())
    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    print("Running:" + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')

    return stdoutdata

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
    PlotDistribution(num_enhancers, "NumberOfGenesPerEnhancer", outdir)
    PlotDistribution(GeneCounts, "NumberOfEnhancersPerGene", outdir)
    PlotDistribution(enhancergeneperchrom, "EnhancersPerChromosome", outdir)
    PlotDistribution(distance, "EnhancerGeneDistance", outdir)

    with open(os.path.join(outdir,"QCSummary.txt"), "w") as f:
        f.write("Median Enhancer Per Gene:")
        f.write(str(GeneMedian))
        f.write("\t")
        f.write(str(GeneStdev))
        f.write("\n")
        f.write("Median Genes Per Enhancer:")
        f.write(str(median_genes_per_enhancer))
        f.write("\t")
        f.write(str(stdev_genes_per_enhancer))
        f.write("\n")
        f.write("Mean Enhancer Per Gene:")
        f.write(str(GeneMean))
        f.write("\n")
        f.write("Mean Genes Per Enhancer:")
        f.write(str(mean_genes_per_enhancer))
        f.write("\n")
        f.write("E-G distance:")
        f.write(str(np.median(distance)))
        f.write("\t")
        f.write(str(np.mean(distance)))
        f.write("\t")
        f.write(str(np.std(distance)))
        f.write("\n")
        f.write("Number of Enhancers/Chrom:")
        f.write(str(median_enhancergeneperchrom))
        f.write("\t")
        f.write(str(mean_enhancergeneperchrom))
        f.write("\t")
        f.write(str(stdev_enhancergeneperchrom))
        f.write("\n")
        f.write("E-G 10th quantile:")
        f.write(str(thquantile))
        f.write("\n")
        f.write("E-G 90th quantile:")
        f.write(str(testthquantile))
        f.close()

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


    

def NeighborhoodFileQC(neighborhood_dir, outdir):

    x = glob.glob(os.path.join(neighborhood_dir, "Enhancers.DHS.*CountReads.bedgraph"))
    y = glob.glob(os.path.join(neighborhood_dir, "Genes.TSS1kb.DHS.*CountReads.bedgraph"))
    z = glob.glob(os.path.join(neighborhood_dir, "Genes.DHS.*CountReads.bedgraph"))
    
    data = pd.read_csv(x[0],sep="\t", header=None)
    data1 = pd.read_csv(y[0],sep="\t", header=None)
    data2 = pd.read_csv(z[0],sep="\t", header=None)

    counts = data.iloc[:, 3].sum()
    counts2 = data1.iloc[:,3].sum()
    counts3 = data2.iloc[:,3].sum()
    with open(os.path.join(outdir, "PeakFileQCSummary.txt"), "a") as f:
        f.write("Counts in Enhancers/GenesTSS/Genes:")
        f.write("\t")
        f.write(str(counts))
        f.write("\t")
        f.write(str(counts2))
        f.write("\t")
        f.write(str(counts3))
        f.write("\n")
        f.close()

# Generates peak file metrics
def PeakFileQC(macs_peaks, outdir):
    if macs_peaks.endswith(".gz"):
        peaks = pd.read_csv(macs_peaks, compression="gzip", sep="\t", header=None)
    else:
        peaks = pd.read_csv(macs_peaks, sep="\t", header=None)

    # Calculate metrics for candidate regions 
    outfile = os.path.join(outdir, os.path.basename(macs_peaks) + ".candidateRegions.bed")
    
    candidateRegions = pd.read_csv(outfile, sep="\t", header=None)
    candidateRegions['dist'] = candidateRegions[2] - candidateRegions[1]
    candreg = list(candidateRegions['dist'])
    PlotDistribution(candreg, 'WidthOfCandidateRegions', outdir)
    
    # Calculate metrics to grab closest TSS to peak 
    annotatedFile = os.path.join(outdir, os.path.basename(macs_peaks) + ".annotated_peaks.bed")
    annotatedPeaks = pd.read_csv(annotatedFile, sep="\t", header=None)
    annote = annotatedPeaks.loc[annotatedPeaks.iloc[:,9]!=-1]
    mean, median, stdev = grabStatistics(annote.iloc[:,9])

#    PlotDistribution(np.array(annotatedPeaks.iloc[:,9]), "DistanceOfPeakToClosestTSS", outdir)

    # Calculate width of peaks 
    peaks['dist'] = peaks[2]-peaks[1]
    peaks_array = list(peaks['dist'])

#    PlotDistribution(peaks_array, "WidthOfPeaks", outdir)
    with open(os.path.join(outdir, "PeakFileQCSummary.txt"),"w") as f:
        f.write(str(macs_peaks))
        f.write("\n")
        f.write("Number of peaks:")
        f.write(str(len(peaks['dist'])))
        f.write("\n")
        f.write("Median width of peak:")
        f.write(str((peaks['dist'].median())))
        f.write("\t")
        f.write(str(peaks['dist'].mean()))
        f.write("\t")
        f.write(str(peaks['dist'].std()))
        f.write("\n")
        f.write("Median Distance of Peak to Closest TSS:")
        f.write(str(median))
        f.write("\t")
        f.write(str(mean))
        f.write("\t")
        f.write(str(stdev))
        f.write("\n")
        f.write("Number of Candidate Regions:")
        f.write(str(len(candidateRegions['dist'])))
        f.write("\n")
        f.write("Median width of Candidate Regions:")
        f.write(str((candidateRegions['dist'].median())))
        f.write("\t")
        f.write(str(candidateRegions['dist'].mean()))
        f.write("\t")
        f.write(str(candidateRegions['dist'].std()))
        f.write("\n")
        f.close()

# Plots and saves a distribution as *.png
def PlotDistribution(array, title, outdir):
    ax = sns.distplot(array)
    ax.set_title(title)
    ax.set_ylabel('Estimated PDF of distribution')
    ax.set_xlabel('Counts')
    fig = ax.get_figure()
    outfile = os.path.join(outdir, str(title)+".pdf")
    fig.savefig(outfile, format='pdf')
