#! bin/python3
import sys
from metrics import *
from tools import *
import glob
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--peaks_outdir", help="File containing peaks directory")
    parser.add_argument("--macs_peaks", help="narrowPeak file output by macs2. Must include summits (--call-summits)")
    parser.add_argument("--neighborhood_outdir", help="File containing neighborhood directory")
    parser.add_argument("--preds_outdir", help="File containing predictions directory")
    parser.add_argument("--quantile", default=False, help="Boolean to plot quantile norm")
    args = parser.parse_args()
    return args 
    

def generateQCMetrics(args):
    genome_tss = "../reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed"
    prediction = "{}/EnhancerPredictionsFull.txt".format(args.preds_outdir)
    grab_nearest_tss_from_peak(args.macs_peaks, genome_tss, args.peaks_outdir)
    # read prediction file
    prediction_df = pd.read_csv(prediction, sep="\t")
    # Generate QC Summary.txt in Predictions Directory
    GrabQCMetrics(prediction_df, preds_outdir)
    # Generate PeakFileQCSummary.txt in Peaks Directory#    
    PeakFileQC(args.macs_peaks, args.peaks_outdir)
    # Appends Percentage Counts in Promoters into PeakFileQCSummary.txt
    NeighborhoodFileQC(args.neighborhood_outdir, args.peaks_outdir)

    if args.quantile:
        # Generate Quantilenorm Plots
        enhancers = "{}/EnhancerList.txt".format(args.neighborhood_outdir)
        EnhancerList = pd.read_csv(enhancers, sep="\t")
        title="_QuantileNorm"
        PlotQuantilePlot(EnhancerList, title, outdir)

if __name__=="__main__":
    args = parse_args()
    generateQCMetrics(args)
