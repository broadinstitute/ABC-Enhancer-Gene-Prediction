#! bin/bash 

while read p;
do
	python /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/data/ENCODEPortal/workflow/scripts/encode_task_filter.py --bam $2/${p} --out-dir $3 $4
#	python encode_task_filter.py --bam $2/$p --out-dir /srv/scratch/kmualim/celltype_files/
done < $1
