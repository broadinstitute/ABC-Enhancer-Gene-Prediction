#! bin/bash 

outdir=$1
exptdata=$2
for num in {0..10};
do
	Rscript comparePredictionsToExperiment.R --predictions $1/train_${num} --experimentalData $2 --plotConfig /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/comparison/CRISPR/plot.config.txt --predConfig /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/comparison/CRISPR/pred.config.txt
	bash renameFiles.sh train_${num} $1	
done	
