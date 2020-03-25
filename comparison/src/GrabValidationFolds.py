import pandas as pd
import numpy as np
import os,sys
import subprocess 
from sklearn import metrics
from scipy import interpolate
from sklearn.model_selection import train_test_split

## sys.argv[1] = name of enhancer-gene prediction file 
## sys.argv[2] = outdir 
## sys.argv[3] = experimental data to compare against 

file_a = sys.argv[1]
outdir = sys.argv[2] 
exptdata = sys.argv[3]

if not os.path.exists(outdir):
        os.makedirs(outdir)
print("Reading Enhancer Predictions File....")
enh = pd.read_csv(file_a, sep="\t", compression="gzip")

print("Splitting Enhancer Predictions file into 10 Cross Validation Fold.....")
for i in range(10):
	print("Splitting for file {}".format(i))
	enh.sample(frac=1)
	train, test = train_test_split(enh, test_size=0.05)
	train.to_csv(os.path.join(outdir, "train_{}".format(i)), sep="\t")

print("Done splitting files!")
print("Running Comparison Script on all 10 Cross Validation Folds ....")
# Run comparison code on all K-fold
rc = subprocess.Popen("./RunComparisonScript.sh {} {}".format(outdir, exptdata), shell=True)

# fill na and return auc
def get_metric(df, recall):
    df_1 = df.fillna(1.0)
    recall = 0.6
    f = interpolate.interp1d(df_1['recall'], df_1['precision'])
    precision = f(recall)
    y = f(df_1['recall'])
    return metrics.auc(df_1['recall'], y), precision

print("Grabbing metric values from pr curves")
auprc_values = []
precision_values = [] 
recall = 0.6
for i in range(10):
	print("Grabbing values from file {}".format(i))
	val = pd.read_csv(os.path.join(outdir, "pr.curve.train_{}.txt".format(i)), sep="\t")
	ABC_val = val.loc[val["pred.col"]=='ABC.Score']
	auprc, precision = get_metric(ABC_val, recall)
	auprc_values.append(auprc)
	precision_values.append(precision)

print("Saving values into final tsv file")
open(os.path.join(outdir,"auprc_precision.tsv"), 'w').write(
            "\n".join([str(i)+"\t"+str(j) for i,j in
                       zip(auprc_values, precision_values) ]))

print("Done!")
# open files 
#command = "bash RunComparisonScript.sh" 
#p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
#print(" \n Running: " + sort_command + "\n")
#(stdoutdata, stderrdata) = p.communicate()
#err = str(stderrdata, 'utf-8')
