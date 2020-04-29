# Running workflow 

This Snakefile directory contains four directories, namely: 

## Existing Directories
envs : contains the software environment in conda + config yaml file for snakemake rules <br>
* working directories need to be configured for the current system <br>
* parameters can all be adjusted in the environment folders <br>
* all the input/output files are very specific to where its located relative to the ABC-rep, so working directories need to be accurately specified 

output : contains the output files from the rules from the metadata snakemake workflow, data input files for input into preprocessing rule + input data lookup into abc 
* most of the output files here are generated from Snakefile in download which runs scripts/grabDownload.py and scripts/download_log.sh : 
	*  ```download_log.sh``` outputs metadata files for all hg19/GHCR38 DHS/ATAC/H3K27ac which is used as input into preprocessing bam files 
	* ```grabDownload.py``` uses the metadata file to perform filtering steps as well as to generate input data lookup table to pair (DHS/H3K27ac) for each biosample term name, and to generate the list of pairedend/singleend bam files for removing duplicates.  
	* Experiments_ToCombine.txt also gets generated and lists the bam files to combine. 	

rules : contains rules for 
1. downloading metadata + bamfiles
2. preprocessing downloaded bam files (handles pairedend + singleend)
3. running abc

* To run each snakefile workflow: simply cd into the directory and run the command ***snakemake***
	* download: download bamfiles 
	* preprocessing: preprocesses bamfiles based on DHS/H3K27ac bamfiles 
	* abc_code : runs ABC on input_data_lookup.txt files and generates enhancer-gene predictions
	
