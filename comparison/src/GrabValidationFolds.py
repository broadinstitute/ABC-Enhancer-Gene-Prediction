import argparse
import pandas as pd
import numpy as np
import os,sys
import subprocess 
from sklearn import metrics
from scipy import interpolate
from sklearn.model_selection import train_test_split

def parse_args():
    parser = argparse.ArgumentParser(description="Script for Cross Validation, to compare ABC Predictions with Experimental data")
    parser.add_argument('--preds', help="Path to EnhancerPredictionsAllPutative.txt.gz file")
    parser.add_argument('--outdir', help="Outdir to save cross validation split files")
    parser.add_argument('--expt', help="Path to Experimental Datafile")
    args = parser.parse_args()
    return args 

def read_file(file_a):
    print("Reading Enhancer Gene Prediction File....")
    enh = pd.read_csv(file_a, sep="\t", compression="gzip")
    return enh

def split_files(enh, outdir):
    print("Splitting Enhancer Gene Predicion File....")
    for i in range(10):
        print("Splitting for file {}".format(i))
        enh.sample(frac=1)
        train, test = train_test_split(enh, test_size=0.05)
        train.to_csv(os.path.join(outdir, "train_{}".format(i)), sep="\t")
    print("Done splitting files!")

def run_comparison_script(outdir, expt_data):
    print("Running comparison script for split files")
    rc = subprocess.call("bash RunComparisonScript.sh {} {}".format(outdir, expt_data), shell=True)

## fill na and return auc
def get_metric(df, recall):
    df_1 = df.fillna(1.0)
    recall = 0.6
    f = interpolate.interp1d(df_1['recall'], df_1['precision'])
    precision = f(recall)
    y = f(df_1['recall'])
    return metrics.auc(df_1['recall'], y), precision

def consolidate_Values(outdir):
    auprc_values = []
    precision_values = [] 
    recall = 0.6
    for i in range(10):
        try:
    	    print("Grabbing values from file {}".format(i))
    	    val = pd.read_csv(os.path.join(outdir, "pr.curve.train_{}.txt".format(i)), sep="\t")
    	    ABC_val = val.loc[val["pred.col"]=='ABC.Score']
    	    auprc, precision = get_metric(ABC_val, recall)
    	    auprc_values.append(auprc)
    	    precision_values.append(precision)
        except:
            continue
    print("Saving values into final tsv file")
    open(os.path.join(outdir,"auprc_precision.tsv"), 'w').write(
                "\n".join([str(i)+"\t"+str(j) for i,j in
                           zip(auprc_values, precision_values) ]))
    
    print("Done!")

if __name__=="__main__":
    args = parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
#    enh = read_file(args.preds)
#    split_files(enh, args.outdir)
    run_comparison_script(args.outdir, args.expt)
    consolidate_Values(args.outdir)

