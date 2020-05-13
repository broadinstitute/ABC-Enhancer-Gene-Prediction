#! bin/bash

####################################################################################
## Code for running ENCODE Snakefiles on Kundaje Lab Cluster 
## Kristy Mualim
## May 11, 2020

CODEDIR=/users/kmualim/updated_ABC/github/; cd $PROJECT
 

##########################################################################
## Run 'download' snakefile
snakemake -s "$CODEDIR/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/rules/download/Snakefile" --configfile $CODEDIR/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/envs/wd.yaml --directory $CODEDIR/ &> logs/download.out

##########################################################################
## Run 'preprocess' snakefile
snakemake -s "$CODEDIR/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/rules/preprocessing/Snakefile" --configfile $CODEDIR/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/envs/wd.yaml --directory $CODEDIR/ &> logs/preprocess.out

##########################################################################
## Run 'abc_code' snakefile
snakemake -s "$CODEDIR/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/rules/abc_code/Snakefile" --configfile $CODEDIR/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/envs/wd.yaml --directory $CODEDIR/ &> logs/abc_code.out
