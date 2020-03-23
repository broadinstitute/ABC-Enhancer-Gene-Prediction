import os,sys
import pandas as pd
import numpy as np
import seaborn as sns
import glob
from subprocess import check_call, check_output, PIPE, Popen, getoutput, CalledProcessError

def grab_nearest_tss_from_peak(macs_peaks, genome_tss, outdir):
<<<<<<< HEAD
    # Grab nearest tss from peak
    outfile = os.path.join(outdir, os.path.basename(macs_peaks))
    files = pd.read_csv(outfile, sep="\t")
    annotated_peaks = os.path.join(outdir, os.path.basename(macs_peaks) + ".annotated_peaks.bed")
    sort_command = "sort -k1,1 -k2,2n {genome_tss} > {genome_tss}.sorted"
    sort_command = sort_command.format(**locals())
    p = Popen(sort_command, stdout=PIPE, stderr=PIPE, shell=True)
    print("Running:" + sort_command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')

    command = "bedtools closest -a {outfile} -b {genome_tss}.sorted -d > {annotated_peaks}"
=======
    
    outfile = os.path.join(outdir, os.path.basename(macs_peaks))
    files = pd.read_csv(outfile, sep="\t")
    annotated_peaks = os.path.join(outdir, os.path.basename(macs_peaks) + ".annotated_peaks.bed")
    
    # command to get closest gene tss using bedtools 
    command = "bedtools closest -a {outfile} -b {genome_tss} -d > {annotated_peaks}"
>>>>>>> 83c5de0ed2e980b872d0eb9e4759ec6eceae265a
    command = command.format(**locals())
    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    print("Running:" + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')

    return stdoutdata

# Generates QC Prediction Metrics:
def GrabQCMetrics(prediction_df, outdir):
    GeneCounts = prediction_df.groupby(['TargetGene']).size()
    GeneCounts.to_csv(os.path.join(outdir,"EnhancerPerGene.txt"), sep="\t")

    # Grab Number of Enhancers Per Gene
    GeneMedian = prediction_df.groupby(['TargetGene']).size().median()
    GeneMean = prediction_df.groupby(['TargetGene']).size().mean()
    GeneStdev = prediction_df.groupby(['TargetGene']).size().std()

    # Grab Number of genes per enhancers
    num_enhancers = prediction_df[['chr', 'start', 'end']].groupby(['chr', 'start', 'end']).size()
    num_enhancers.to_csv(os.path.join(outdir,"GenesPerEnhancer.txt"), sep="\t")
    mean_genes_per_enhancer = prediction_df[['chr', 'start', 'end']].groupby(['chr', 'start', 'end']).size().mean()
    stdev_genes_per_enhancer = prediction_df[['chr', 'start', 'end']].groupby(['chr', 'start', 'end']).size().std()
    median_genes_per_enhancer = prediction_df[['chr', 'start', 'end']].groupby(['chr', 'start', 'end']).size().median()
    
    # Grab Number of Enhancer-Gene Pairs Per Chromsome
    enhancergeneperchrom = prediction_df.groupby(['chr']).size()
    enhancergeneperchrom.to_csv(os.path.join(outdir, "EnhancerGenePairsPerChrom.txt"), sep="\t")
    mean_enhancergeneperchrom = prediction_df.groupby(['chr']).size().mean()
    stdev_enhancergeneperchrom = prediction_df.groupby(['chr']).size().std()
    median_enhancergeneperchrom = prediction_df.groupby(['chr']).size().median()

    # Enhancer-Gene Distancee
    distance = np.array(prediction_df['distance'])
    thquantile = np.percentile(distance, 10)
    testthquantile = np.percentile(distance, 90)
<<<<<<< HEAD

    # Plot Distributions and save as png
=======
    
    # Plot Distributions and save as pdf
>>>>>>> 83c5de0ed2e980b872d0eb9e4759ec6eceae265a
    PlotDistribution(num_enhancers, "NumberOfGenesPerEnhancer", outdir)
    PlotDistribution(GeneCounts, "NumberOfEnhancersPerGene", outdir)
    PlotDistribution(enhancergeneperchrom, "EnhancersPerChromosome", outdir)
    PlotDistribution(distance, "EnhancerGeneDistance", outdir)

    # save file in predictions directory 
    with open(os.path.join(outdir,"QCSummary.txt"), "w") as f:
        f.write("Enhancer Per Gene:")
        f.write(str(GeneMedian))
        f.write("\t")
        f.write(str(GeneMean))
        f.write("\t")
        f.write(str(GeneStdev))
        f.write("\n")
        f.write("Genes Per Enhancer:")
        f.write(str(median_genes_per_enhancer))
        f.write("\t")
        f.write(str(mean_genes_per_enhancer))
        f.write("\t")
        f.write(str(stdev_genes_per_enhancer))
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
<<<<<<< HEAD


=======
    
# Grab Number of Counts that lie in Enhancers, GeneTSS and Gene Bodies - appended to PeakFileQCSummary.txt 
>>>>>>> 83c5de0ed2e980b872d0eb9e4759ec6eceae265a
def NeighborhoodFileQC(neighborhood_dir, outdir):
    x = glob.glob(os.path.join(neighborhood_dir, "Enhancers.DHS.*"))
    y = glob.glob(os.path.join(neighborhood_dir, "Genes.TSS1kb.DHS.*"))
    z = glob.glob(os.path.join(neighborhood_dir, "Genes.DHS.*"))
    
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
        f.close()

# Generates peak file metrics
def PeakFileQC(macs_peaks, outdir):
    if macs_peaks.endswith(".gz"):
        peaks = pd.read_csv(macs_peaks, compression="gzip", sep="\t", header=None)
    else:
        peaks = pd.read_csv(macs_peaks, sep="\t", header=None)
    outfile = os.path.join(outdir, os.path.basename(macs_peaks) + ".candidateRegions.bed")
    
    # Grab width of candidate regions
    candidateRegions = pd.read_csv(outfile, sep="\t", header=None)
    candidateRegions['dist'] = candidateRegions[2] - candidateRegions[1]
    candreg = list(candidateRegions['dist'])
    PlotDistribution(candreg, 'WidthOfCandidateRegions', outdir)

    # Grab Distance of Peak To Closest TSS
    annotatedFile = os.path.join(outdir, os.path.basename(macs_peaks) + ".annotated_peaks.bed")
    annotatedPeaks = pd.read_csv(annotatedFile, sep="\t", header=None)
    median = annotatedPeaks.iloc[:, 9].median()
    mean = annotatedPeaks.iloc[:, 9].mean()
    stdev = np.std(np.array(annotatedPeaks.iloc[:, 9]))
    PlotDistribution(np.array(annotatedPeaks.iloc[:,9]), "DistanceOfPeakToClosestTSS", outdir)

    # Grab Width of Peaks
    peaks['dist'] = peaks[2]-peaks[1]
    peaks_array = list(peaks['dist'])
    PlotDistribution(peaks_array, "WidthOfPeaks", outdir)

    # Save to PeakFileQCSummary
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
        f.write(str((peaks['dist'].median())))
        f.write("\t")
        f.write(str(peaks['dist'].mean()))
        f.write("\t")
        f.write(str(peaks['dist'].std()))
        f.write("\n")
        f.close()

# Plots and saves a distribution as *.pdf
def PlotDistribution(array, title, outdir):
    ax = sns.distplot(array, hist=True)
    ax.set_title(title)
    ax.set_ylabel('Estimated PDF of distribution')
    ax.set_xlabel(str(title))
    fig = ax.get_figure()
    outfile = os.path.join(outdir, str(title)+".pdf")
    fig.savefig(outfile, format='pdf')
    plt.clf()
