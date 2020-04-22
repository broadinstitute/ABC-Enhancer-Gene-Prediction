This Snakefile directory contains four directories, namely: 

envs : contains the software environment in conda + config yaml file for snakemake rules
	- working directories need to be configured for the current system 
	- parameters can all be adjusted in the environment folders 
	- all the input/output files are very specific to where its located relative to the ABC-rep, so working directories need to be accurately specified 

output : contains the output files from the rules from the metadata snakemake workflow, data input files for input into preprocessing rule + input data lookup into abc 

rules : contains rules for downloading metadata + bamfiles, preprocessing downloaded bam files (handles pairedend + singleend), running abc
	- To run each snakefile workflow: simply cd into the directory and run the command """ snakemake """
