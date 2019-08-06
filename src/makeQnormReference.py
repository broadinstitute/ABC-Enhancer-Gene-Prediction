#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Run neighborhood for a given cell type',
                                     epilog=epilog,
                                     formatter_class=formatter)

    parser.add_argument('--enhancers', required=required_args, help="EnhancerList file")
    parser.add_argument('--outfile', default="", help="columns for which to compute qnorm reference")
    parser.add_argument('--cols', default="DHS.RPM,H3K27ac.RPM", help="columns for which to compute qnorm reference")

    args = parser.parse_args()
    return args


def makeQnorm(args):
    enhancers = pd.read_csv(args.enhancers, sep="\t")

    ref = []

    # import pdb
    # pdb.set_trace()

    for enh_type in ['any','promoter','nonpromoter']:
        if enh_type == 'any':
            this_enhancers = enhancers
        elif enh_type == 'promoter':
            this_enhancers = enhancers.loc[np.logical_or(enhancers['class'] == "tss", enhancers['class'] == "promoter")]
        elif enh_type == 'nonpromoter':
            this_enhancers = enhancers.loc[np.logical_and(enhancers['class'] != "tss" , enhancers['class'] != "promoter")]
        else:
            error("Wrong type")

        this_ref = pd.DataFrame({'enh_class' : enh_type, 'quantile' : np.concatenate([np.linspace(.01,.99,99), np.linspace(.991,.999,9), np.linspace(.9991,.9999,9)])})
        this_ref['rank'] = ((1 - this_ref['quantile']) * this_enhancers.shape[0]).round(decimals=0)

        cols = set(set(vars(args)['cols'].split(",")) & set(this_enhancers.columns))
        for col in cols:
            this_ref[col] = np.percentile(this_enhancers[col].values, this_ref['quantile'].values*100)

        ref.append(this_ref)

    # ref = pd.DataFrame({'quantile' : np.concatenate([np.linspace(0,.99,100), np.linspace(.991,.999,9), np.linspace(.9991,.9999,9)])})
    # ref['rank'] = ((1 - ref['quantile']) * enhancers.shape[0]).round(decimals=0)

    # cols = set(set(vars(args)['cols'].split(",")) & set(enhancers.columns))
    # for col in cols:
    #     ref[col] = np.percentile(enhancers[col].values, ref['quantile'].values*100)

    ref = pd.concat(ref)
    ref.to_csv(args.outfile, sep="\t", index=False, float_format="%.5f")

def main(args):
    makeQnorm(args)

if __name__ == '__main__':
    args = parseargs()
    main(args)
