import yaml
ASSAYS=["DHS", "H3K27ac", "ATAC"]
endtype = ["pairend", "singleend"]

config = yaml.load(open('../envs/wd.yaml'))
genome = config['params_download']['genome']

rule all:
	input: expand("bam_{genome_build}_{assay}.tsv", assay=ASSAYS, genome_build=genome), expand("{pairtype}_bam_{genome_build}_{assays}.tsv", assays=ASSAYS, pairtype=endtype, genome_build=genome)

rule download_bam:
	input:
		exe = config['prefix']['workflow_scripts']+"log.sh"
        message: "Executing log.sh to download files"
	output: 
		outdir = config['params_preprocess']['preprocess_dir']
		files = expand("bam_{genome_build}_{assay}.tsv", assay=ASSAYS, genome_build=genome)
	log: "logs/download.log"
	shell: "bash {input.exe} {} {output.outdir}".format(genome)

rule obtain_paired_single_files:
	input:
		exe = config['prefix']['workflow_scripts']+"log_kristy.sh"
	output: 
		outdir = config['params_preprocess']['preprocess_dir']
		files =expand("{pairtype}_bam_{genome_build}_{assays}.tsv", assays=ASSAYS, pairtype=endtype, genome_build=genome)
	message: "Executing log_kristy.sh to download paired-end / single-end file information"
	log: "logs/pe_se_download.log"
	shell: "bash {input.exe} {} {output.outdir}".format(genome)
