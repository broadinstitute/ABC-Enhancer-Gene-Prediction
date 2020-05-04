#! bin/python3

import sys, os
import pandas as pd
import numpy

input_file=sys.argv[1]
pvalue=sys.argv[2]
outdir=sys.argv[3]
input_dir=sys.argv[4]


if __name__ == '__main__':
    args = parse_arguments()
    input_files = pd.read_csv(args.input_file, sep="\t", header=None)
    if args.apply_threading:
        threads = args.threads
        with Pool(int(threads)) as p:
            p.map(main, input_files[0])
    else:

while read p;
def main():
    a=($p)
    echo ${a[0]} ${a[1]}
	macs2 callpeak -f BAM -g hs -p $pvalue --call-summits --outdir $outdir/Peaks_${a[0]} -t $input_dir/${a[1]}
done < $input_file

