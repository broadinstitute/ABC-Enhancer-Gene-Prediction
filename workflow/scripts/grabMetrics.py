#! bin/python3
import sys
import pandas as pd
from metrics import *
from tools import *
import glob
import argparse
import pickle

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--macs_peaks", help="narrowPeak file output by macs2. (eg. /users/kmualim/K562/macs2_peaks.narrowPeak)")
    parser.add_argument("--preds_file", help="Prediction file output by ABC")
    parser.add_argument("--hicDir", help="HiC Directory")
    parser.add_argument('--hic_powerlaw_fit', help="hic.powerlaw.csv file generated from compute_powerlaw_fit_from_hic.py")
    parser.add_argument("--neighborhood_outdir", help="Neighborhood Directory")
    parser.add_argument("--outdir", help="Predictions Directory")
    args = parser.parse_args()
    return args 
    

def generateQCMetrics(args):
    #prediction = "{}/EnhancerPredictionsFull.txt".format(args.preds_outdir)
   # read prediction file
    prediction_df = pd.read_csv(args.preds_file, sep="\t")
   # Generate QC Summary.txt in Predictions Directory
    pred_metrics = GrabQCMetrics(prediction_df, args.outdir)
   # Generate PeakFileQCSummary.txt in Peaks Directory#    
    pred_metrics = PeakFileQC(pred_metrics, args.macs_peaks, args.outdir)
    # Appends Percentage Counts in Promoters into PeakFileQCSummary.txt
    pred_metrics = NeighborhoodFileQC(pred_metrics, args.neighborhood_outdir, args.outdir, "DHS")
    pred_metrics = NeighborhoodFileQC(pred_metrics, args.neighborhood_outdir, args.outdir, "H3K27ac")
    # Generate Powerlaw QC plot from HiC directory
    if args.hicDir:
        assert args.hic_powerlaw_fit, "Must provide powerlaw fit with HiC data"
        powerlaw_fit = pd.read_csv(args.hic_powerlaw_fit, sep="\t")
        HiCQC(prediction_df, powerlaw_fit, args.outdir)
    
    with open("{}/QCSummary.p".format(args.outdir),"wb") as f:
        pickle.dump(pred_metrics, f) 
    

if __name__=="__main__":
    args = parse_args()
    generateQCMetrics(args)
