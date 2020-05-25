####################################################################################
## Code for running ENCODE Snakefiles on Broad Server (gsa queue)
## Jesse Engreitz
## May 11, 2020

PROJECT=/seq/lincRNA/RAP/ABC/200330_ENCODE/200525; cd $PROJECT
CODEDIR=$PROJECT/ABC-Enhancer-Gene-Prediction/

#Activate conda environment
conda activate /seq/lincRNA/RAP/ABC/200330_ENCODE/env/abcenv

#source /broad/software/scripts/useuse
# reuse .python-3.5.1 .samtools-1.8 .tabix-0.2.5 BEDTools
# source /seq/lincRNA/Ben/VENV_MIP/bin/activate

#SAMPLESHEET=$1
#CONFIG=$2 
# example config: /seq/lincRNA/Projects/MIP/mip-pipeline/file_input/config.json

mkdir -p logs download 


##########################################################################
## Run 'download' snakefile
snakemake -s "$CODEDIR/Snakefiles/workflow/rules/download/Snakefile" \
	--configfile $CODEDIR/Snakefiles/workflow/envs/wd.yaml \
	--config working_dir=$PROJECT download_dir=$PROJECT/download output_data_dir=$PROJECT/download \
	 &> logs/download.out &
	#--profile $CODEDIR/Snakefiles/workflow/envs/BroadGSA/ \  ## Don't include cluster job submission profile -- want to run this pipeline locally


##########################################################################
## Run 'preprocess' snakefile
## TO EDIT Joe
snakemake -s "$CODEDIR/Snakefiles/workflow/rules/preprocess/Snakefile" \
	--profile $CODEDIR/Snakefiles/workflow/envs/BroadGSA/ \
	--configfile $CODEDIR/Snakefiles/workflow/envs/wd.yaml \
	## EDIT THIS: --config working_dir=$PROJECT download_dir=$PROJECT/download output_data_dir=$PROJECT/download \
	 &> logs/download.out &


##########################################################################
## Run 'ABC' snakefile
## TO EDIT Joe
snakemake -s "$CODEDIR/Snakefiles/workflow/rules/abc_code/Snakefile" \
	--configfile $CODEDIR/Snakefiles/workflow/envs/wd_Broad.yaml \
	--config working_dir=$PROJECT download_dir=$PROJECT/download output_data_dir=$PROJECT/download \
			predictions_results_dir=$PROJECT/abc_out \
			hic_dir=/seq/lincRNA/RAP/External/Rao2014-HiC/190924_new_structure/hic/191020_Average/out/ \
		#--profile $CODEDIR/Snakefiles/workflow/envs/BroadGSA/ \
	 &> logs/download.out &
