#! bin/bash 

input_file=$1
chrom_sizes=$2
datadir=$3
codedir=$4
outdir=$5
genes=$6
ubiquitous_genes=$7
qnorm=$8

while read p;
do
	a=($p)
	echo ${a[0]} ${a[1]}
	python $codedir/src/run.neighborhoods.py --candidate_enhancer_regions $outdir/Peaks_${a[0]}/NA_peaks.narrowPeak.sorted.candidateRegions.bed --DHS $datadir/${a[1]} --H3K27ac $datadir/${a[2]} --genes $genes --ubiquitously_expressed_genes $ubiquitous_genes --cellType ${a[0]} --outdir $outdir/Neighborhoods_${a[0]} --chrom_sizes $chrom_sizes --qnorm $qnorm

done < $input_file

