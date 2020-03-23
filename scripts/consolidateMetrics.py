import os, sys
import pandas as pd
import numpy as np

celltypes = sys.argv[1]
celltypes_files = pd.read_csv(celltypes, sep="\t", header=None)

def fill_peakfile(peaks):
    title_peaks_metric = [str(peaks_metric[peak][0]).split(":")[0] for peak in range(len(peaks_metric)) if peak in range(1,6)]
    total_peaks_metric = [str(peaks_metric[peak][0]).split(":")[1] for peak in range(len(peaks_metric)) if peak in range(1,6)]
    percentage_counts = peaks_metric[-1]
    total_sum = np.sum([int(percentage_counts[i]) for i in range(1,4)])
    title_peaks_metric.append("PercentageCountsInPromoterRegions")
    total_peaks_metric.append(int(percentage_counts[2])/total_sum)
    return title_peaks_metric, total_peaks_metric

def fill_predsfile(preds_metric):
    title_preds_metric = [str(preds_metric[peak][0]).split(":")[0] for peak in range(len(preds_metric)) if peak in range(0,6)]
    total_preds_metric = [str(preds_metric[peak][0]).split(":")[1] for peak in range(len(preds_metric)) if peak in range(0,6)]
    return title_preds_metric, total_preds_metric
    
def update_dict(peak_dict_value, value, dictionary, title):
    peak_dict_values.append(peaks)
    total_peaks[title] = peak_dict_values
    return total_peaks
    

def aggregate(celltypes_files):
    total_peaks_metric_dict ={}
    total_preds_metric_dict = {}
    for cell in celltypes_files[0]:
            print(cell)
    peaks_metric = []
    preds_metric = []
    with open("/mnt/lab_data3/kmualim/PeakAndNeighborhoods/Peaks_{}/PeakFileQCSummary.txt".format(str(cell)), "r") as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            peaks_metric.append(row)
    title_peaks_metric, total_peaks_metric = fill_peakfile(peaks_metric)
    
    with open("/mnt/lab_data3/kmualim/PredictionFiles_qnorm/Predictions_{}_qnorm/QCSummary.txt".format(str(cell)), "r") as file:
        file_reader = csv.reader(file, delimiter='\t')
        for row in file_reader:
            preds_metric.append(row)
    title_preds_metric, total_preds_metric = fill_predsfile(preds_metric)
    
    for title, peaks, title_preds, preds in zip(title_peaks_metric, total_peaks_metric, title_preds_metric, total_preds_metric):
            if cell==celltypes[0][1]:
                val=[]
                val.append(total_peaks_metric_dict[title])
                val.append(peaks)
                total_peaks_metric_dict[title] = val
                
                val=[]
                val.append(total_preds_metric_dict[title_preds])
                val.append(preds)
                total_preds_metric_dict[title_preds] = val
            
            elif cell==celltypes[0][0]:
                total_peaks_metric_dict[title] = peaks
                total_preds_metric_dict[title_preds] = preds
            else:
                val = [i for i in total_peaks_metric_dict[title]]
                val.append(peaks)
                total_peaks_metric_dict[title] = val
                
                val = [i for i in total_preds_metric_dict[title_preds]]
                val.append(preds)
                total_preds_metric_dict[title_preds] = val
     return total_preds_metric_dict, total_peaks_metric_dict 

def merge_dataframes(total_preds_metric_dict, total_peaks_metric_dict, celltype_files, outfile):
    df2 = df.from_dict(total_peaks_metric_dict)
    df2['CellTypes'] = celltype_files[0]

    df = df.from_dict(total_preds_metric_dict)
    df['CellTypes'] = celltype_files[0]
    df3 = df.merge(df2, right_index=True, left_index=True)
    df3.to_csv(outfile, sep="\t", index=False)


def main():
    outfile = "Overall_QCMetrics.txt"
    total_preds_metric_dict, total_peaks_metric_dict = aggregate(celltypes_files)
    merge_dataframes(total_preds_metric_dict, total_peaks_metric_dict, celltype_files, outfile)
    