import pandas as pd
import numpy as np
import os, sys
import pickle

# ../../../out/finalized/DNase_H3K27ac_Kristy_Jill.tsv
# ../../../out/finalized_metadata/DNase_Kristy_Jill.tsv

# example command: python concatenate_metrics.py biosamples.tsv output.tsv

# sys.argv[1] takes in your list of biosamples that you want to concatenate
features = pd.read_csv(sys.argv[1], sep="\t")
biosamples = features["biosample"]
# have the outdir point to where all your predictions file lie
outdir = sys.argv[3]
concatenated = None
for biosample in biosamples:
    print(biosample)
    if os.path.exists("{}/{}/Predictions/QCSummary.p".format(outdir, biosample)):
        filename = pickle.load(
            open("{}/{}/Predictions/QCSummary.p".format(outdir, biosample), "rb")
        )
        data = pd.DataFrame.from_dict(filename, orient="index")
        dataT = data.T
        dataT["biosample"] = biosample
        if concatenated is None:
            concatenated = dataT
        else:
            concatenated = pd.concat([concatenated, dataT])

concatenated.to_csv(sys.argv[2], sep="\t", index=False)
