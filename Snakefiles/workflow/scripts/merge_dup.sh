#! bin/bash 

workingdir="/users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/data/ENCODEPortal/workflow/scripts"
indir=$1
outdir=$2
input_file=$3
threads=$4

while read p;
do
	a=($p)
	echo ${a[1]} ${a[2]} ${a[3]}
	BASENAME="${a[1]%%.*}"
	if [ ${#a[@]} -gt 3 ]
	then
		samtools merge $outdir/${BASENAME}_pooled.nodup.bam $indir/${a[1]} $indir/${a[2]} $indir/${a[3]} --threads $threads
	else
		samtools merge $outdir/${BASENAME}_pooled.nodup.bam $indir/${a[1]} $indir/${a[2]} --threads $threads
	fi
done < 	$input_file

