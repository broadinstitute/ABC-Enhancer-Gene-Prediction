import glob
import os
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


def grabStatistics(arr_values):
    mean = arr_values.mean()
    median = arr_values.median()
    std = arr_values.std()
    return mean, median, std


def sort_by_chrom_order(series: pd.Series, chrom_order: List[str]) -> pd.Series:
    new_series = pd.Series()
    for chr in chrom_order:
        if chr in series:
            new_series[chr] = series[chr]
    return new_series


# Generates QC Prediction Metrics:
def GrabQCMetrics(prediction_df, chrom_order, outdir, pdf_writer):
    EnhancerPerGene = prediction_df.groupby(["TargetGene"]).size()
    EnhancerPerGene.to_csv(os.path.join(outdir, "EnhancerPerGene.tsv"), sep="\t")
    pdf_writer.savefig(
        PlotDistribution(
            EnhancerPerGene,
            "Number of Enhancers Per Gene",
            x_label="Enhancers per Gene",
            density_line=False,
        )
    )
    # Grab Number of Enhancers Per Gene
    GeneMean, GeneMedian, GeneStdev = grabStatistics(EnhancerPerGene)

    # Grab Number of genes per enhancers
    NumGenesPerEnhancer = (
        prediction_df[["chr", "start", "end"]].groupby(["chr", "start", "end"]).size()
    )
    NumGenesPerEnhancer.to_csv(os.path.join(outdir, "GenesPerEnhancer.tsv"), sep="\t")
    pdf_writer.savefig(
        PlotDistribution(
            NumGenesPerEnhancer,
            "Number Of Genes Per Enhancer",
            x_label="Genes Per Enhancer",
            density_line=False,
        )
    )
    (
        mean_genes_per_enhancer,
        median_genes_per_enhancer,
        stdev_genes_per_enhancer,
    ) = grabStatistics(NumGenesPerEnhancer)

    # Grab Number of Enhancer-Gene Pairs Per Chromsome
    enhancergeneperchrom = prediction_df.groupby(["chr"]).size()
    enhancergeneperchrom = sort_by_chrom_order(enhancergeneperchrom, chrom_order)
    enhancergeneperchrom.to_csv(
        os.path.join(outdir, "EnhancerGenePairsPerChrom.txt"), sep="\t"
    )
    pdf_writer.savefig(
        plotBarPlot(
            enhancergeneperchrom,
            enhancergeneperchrom.index,
            "Enhancer-Gene Pairs Per Chromosome",
            x_label="Number of E-G pairs",
            y_label="Chromosome",
        )
    )
    (
        mean_enhancergeneperchrom,
        median_enhancergeneperchrom,
        stdev_enhancergeneperchrom,
    ) = grabStatistics(enhancergeneperchrom)

    # Enhancer-Gene Distancee
    distance = prediction_df["distance"]
    distance = distance[distance > 0]
    log_dist = distance.apply(np.log10)
    pdf_writer.savefig(
        PlotDistribution(
            log_dist,
            "Enhancer-Gene Distance",
            x_label="Log10 Distance",
            stat="density",
        )
    )
    thquantile = np.percentile(distance, 10)
    testthquantile = np.percentile(distance, 90)

    # Number of Enhancers
    numEnhancers = len(prediction_df[["chr", "start", "end"]].drop_duplicates())
    pred_metrics = {}
    pred_metrics["MedianEnhPerGene"] = GeneMedian
    pred_metrics["StdEnhPerGene"] = GeneStdev
    pred_metrics["MedianGenePerEnh"] = median_genes_per_enhancer
    pred_metrics["StdGenePerEnh"] = stdev_genes_per_enhancer
    pred_metrics["MeanEnhPerGene"] = GeneMean
    pred_metrics["MeanGenePerEnh"] = mean_genes_per_enhancer
    pred_metrics["MedianEGDist"] = np.median(distance)
    pred_metrics["MeanEGDist"] = np.mean(distance)
    pred_metrics["StdEGDist"] = np.std(distance)
    pred_metrics["NumEnhancersPerChrom"] = median_enhancergeneperchrom
    pred_metrics["MeanEnhancersPerChrom"] = mean_enhancergeneperchrom
    pred_metrics["NumEnhancers"] = numEnhancers
    pred_metrics["EG10th"] = thquantile
    pred_metrics["EG90th"] = testthquantile
    return pred_metrics


def plotBarPlot(x_data, y_data, title, x_label, y_label, color="blue"):
    plt.clf()  # make sure we're starting with a fresh plot
    ax = sns.barplot(x=x_data, y=y_data, color=color, orient="h")
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    return ax.get_figure()


def NeighborhoodFileQC(pred_metrics, neighborhood_dir, feature):
    x = glob.glob(
        os.path.join(
            neighborhood_dir, "Enhancers.{}.*CountReads.bedgraph".format(feature)
        )
    )
    y = glob.glob(
        os.path.join(
            neighborhood_dir, "Genes.TSS1kb.{}.*CountReads.bedgraph".format(feature)
        )
    )
    z = glob.glob(
        os.path.join(neighborhood_dir, "Genes.{}.*CountReads.bedgraph".format(feature))
    )

    data = pd.read_csv(x[0], sep="\t", header=None)
    data1 = pd.read_csv(y[0], sep="\t", header=None)
    data2 = pd.read_csv(z[0], sep="\t", header=None)

    counts = data.iloc[:, 3].sum()
    counts2 = data1.iloc[:, 3].sum()
    counts3 = data2.iloc[:, 3].sum()

    pred_metrics["countsEnhancers_{}".format(feature)] = counts
    pred_metrics["countsGeneTSS_{}".format(feature)] = counts2
    pred_metrics["countsGenes_{}".format(feature)] = counts3
    return pred_metrics


# Generates peak file metrics
def PeakFileQC(pred_metrics, macs_peaks, pdf_writer):
    if macs_peaks.endswith(".gz"):
        peaks = pd.read_csv(macs_peaks, compression="gzip", sep="\t", header=None)
    else:
        peaks = pd.read_csv(macs_peaks, sep="\t", header=None)

    # Calculate metrics for candidate regions
    candidateRegions = pd.read_csv(macs_peaks, sep="\t", header=None)
    candidateRegions["dist"] = candidateRegions[2] - candidateRegions[1]
    candreg = candidateRegions["dist"]
    pdf_writer.savefig(
        PlotDistribution(candreg, "WidthOfCandidateRegions", x_label="Width")
    )

    # Calculate width of peaks
    peaks["dist"] = peaks[2] - peaks[1]
    pred_metrics["NumPeaks"] = len(peaks["dist"])
    pred_metrics["MedWidth"] = peaks["dist"].median()
    pred_metrics["MeanWidth"] = peaks["dist"].mean()
    pred_metrics["StdWidth"] = peaks["dist"].std()
    pred_metrics["NumCandidate"] = len(candidateRegions["dist"])
    pred_metrics["MedWidthCandidate"] = candidateRegions["dist"].median()
    pred_metrics["MeanWidthCandidate"] = candidateRegions["dist"].mean()
    pred_metrics["StdWidthCandidate"] = candidateRegions["dist"].std()

    return pred_metrics


# Plots and saves a distribution as *.png
def PlotDistribution(
    pd_series: pd.Series, title, x_label, stat="count", density_line=True
):
    plt.clf()  # make sure we're starting with a fresh plot
    ax = sns.histplot(
        pd_series, kde=density_line, bins=50, kde_kws=dict(cut=3), stat=stat
    )
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(stat.capitalize())

    mean, median, _ = grabStatistics(pd_series)
    mean = round(mean, 3)  # 3 decimal places
    median = round(median, 3)  # 3 decimal places
    plt.axvline(x=mean, color="red", linestyle="-", label=f"Mean={mean}")
    plt.axvline(x=median, color="green", linestyle="-", label=f"Median={median}")
    plt.legend()

    return ax.get_figure()


def HiCQC(df, gamma, scale, pdf_writer):
    # filter for e-g distances of >10kb and <1Mb
    df = df.loc[(df["distance"] > 10000) & (df["distance"] < 1000000)]

    max_samples = 10000
    df = df.sample(min(max_samples, len(df)))
    if len(df):
        pdf_writer.savefig(
            PlotPowerLawRelationship(
                df, "distance", "hic_contact", "E-G Pair HiC Powerlaw Fit", gamma, scale
            )
        )


def PlotPowerLawRelationship(df, x_axis_col, y_axis_col, title, gamma, scale):
    plt.clf()  # make sure we're starting with a fresh plot
    # filter out zeros
    df = df[df[x_axis_col] > 0]
    df = df[df[y_axis_col] > 0]

    log_x_axis_label = f"natural log ({x_axis_col})"
    log_y_axis_label = f"natural log ({y_axis_col})"
    log_x_vals = np.log(df[x_axis_col])
    log_y_vals = np.log(df[y_axis_col])
    values = np.vstack([log_x_vals, log_y_vals])
    kernel = stats.gaussian_kde(values)(values)
    ax = sns.scatterplot(x=log_x_vals, y=log_y_vals, c=kernel, cmap="viridis")

    fitted_y_vals = scale + -1 * gamma * log_x_vals
    sns.lineplot(
        x=log_x_vals, y=fitted_y_vals, color="red", label="Fitted Powerlaw Fit"
    )
    ax.set(title=title)
    ax.set_xlabel(log_x_axis_label)
    ax.set_ylabel(log_y_axis_label)
    return ax.get_figure()
