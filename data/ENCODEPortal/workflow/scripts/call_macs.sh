#! bin/bash 

input_file=$1
pvalue=$2
outdir=$3
input_dir=$4

while read p;
do
	a=($p)
	echo ${a[0]} ${a[1]}
	macs2 callpeak -f BAM -g hs -p $pvalue --call-summits --outdir $outdir/Peaks_${a[0]} -t $input_dir/${a[1]}
done < $input_file

