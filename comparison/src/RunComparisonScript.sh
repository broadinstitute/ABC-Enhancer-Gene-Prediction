#! bin/bash 

outdir = $1
exptdata = $2
for num in {1..10};
do
	Rscript comparePredictionsToExperiment.R --predictions $1/train_${num} --experimentalData $2 --plotConfig CRISPR/plot.config.txt --predConfig CRISPR/pred.config.txt
	bash renameFiles.sh train_${num} $1	
done	
