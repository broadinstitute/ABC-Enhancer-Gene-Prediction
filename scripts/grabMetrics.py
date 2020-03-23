#! bin/python3
import sys
from metrics import *
from tools import *

if __name__=="__main__":
    files = sys.argv[1]
    cells = pd.read_csv(files, sep="\t", header=None)
    
    for i in celltypes: 
        print(i)
        macs_peaks = "/mnt/lab_data3/kmualim/PeakAndNeighborhoods/Peaks_{}/NA_peaks.narrowPeak".format(str(i))#[:-6]))
        genome_tss = "/users/kmualim/updated_ABC/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed"
        peak_outdir = "/mnt/lab_data3/kmualim//PeakAndNeighborhoods/Peaks_{}".format(str(i))#[:-6]))
        outdir="/mnt/lab_data3/kmualim/PredictionFiles_qnorm/Predictions_{}_qnorm".format(i)
        neighborhood_dir = "/mnt/lab_data3/kmualim/NeighborhoodFiles_qnorm/Neighborhoods_{}_qnorm/".format(i)
        prediction="/mnt/lab_data3/kmualim/PredictionFiles_qnorm/Predictions_{}_qnorm/EnhancerPredictionsFull.txt".format(i)
        
        grab_nearest_tss_from_peak(macs_peaks, genome_tss, peak_outdir)
        prediction_df = pd.read_csv(prediction, sep="\t")
        GrabQCMetrics(prediction_df, outdir)
        peak_outdir="/srv/scratch/kmualim/ABC_data/DiffPeakFiles/pval_0_005/Peaks_Predictions_{}_005".format(i)
        PeakFileQC(macs_peaks, peak_outdir)
        NeighborhoodFileQC(neighborhood_dir, peak_outdir)
        enhancers = "/mnt/lab_data3/kmualim/NeighborhoodFiles_qnorm/Neighborhoods_{}/EnhancerList.txt".format(i)
        EnhancerList = pd.read_csv(enhancers, sep="\t")
        title="_QuantileNorm"
        PlotQuantilePlot(EnhancerList, title, outdir)
