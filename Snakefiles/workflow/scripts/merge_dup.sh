#! bin/bash 

workingdir="/users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/data/ENCODEPortal/workflow/scripts"
indir=$1
outdir=$2
input_file=$3
while read p;
do
	a=($p)
	echo ${a[1]} ${a[2]} ${a[3]}
	samtools merge $outdir/${a[1]}_pooled_nodup.bam $indir/${a[1]} $indir/${a[2]} $indir/${a[3]}
	#java -jar AddCommentsToBam.jar INPUT=$outdir/${a[0]}_pooled_nodup.bam.tmp OUTPUT=$outdir/${a[0]}_pooled_nodup.bam	COMMENT=\"ID:samtools	VN:1.2.1	CL: "$1" COMMENT="Picard AddCommentsToBam.jar command was used to add this @CO line and make bam index."	CREATE_INDEX=true CREATE_MD5_FILE=true VALIDATION_STRINGENCY=LENIENT
done < 	$input_file

