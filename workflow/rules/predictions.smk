def _get_run_predictions_hic_params(wildcards):
	hic_dir = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_dir"]
	hic_type = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_type"]
	hic_resolution = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_resolution"]
	if hic_dir:
		return f"--hic_dir {hic_dir} --hic_type {hic_type} --hic_resolution {hic_resolution}"
	else:
		return "--score_column powerlaw.Score"

### run predictions: takes in EnhancerList.txt and GeneList.txt generated from rule call_neighborhoods above and generates Enhancer-Gene Predictions and links
rule run_predictions:
	input:
		enhancers = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
		genes = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
		powerlaw_params_tsv = get_hic_powerlaw_fit_file
	params:
		cellType = lambda wildcards: wildcards.biosample, 
		output_dir = lambda wildcards: os.path.join(RESULTS_DIR, wildcards.biosample, "Predictions"),
		hic_params = _get_run_predictions_hic_params,
		chrom_sizes = config['chrom_sizes'],
		threshold = config['params_predict']['threshold'],
		flags = config['params_predict']['flags'],
		accessibility_feature = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'default_accessibility_feature'],
	conda:
		"../envs/abcenv.yml"
	output: 
		allPutative = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.txt.gz"),
		enhPredictions = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictions.tsv"),
		enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsFull.tsv"),
	resources:
		mem_gb=128,
		runtime_hr=6
	shell:
		"""
		python workflow/scripts/predict.py \
			--enhancers {input.enhancers} \
			--outdir {params.output_dir} \
			{params.hic_params} \
			--powerlaw_params_tsv {input.powerlaw_params_tsv} \
			--chrom_sizes {params.chrom_sizes} \
			--accessibility_feature {params.accessibility_feature} \
			--threshold {params.threshold} \
			--cellType {params.cellType} \
			--genes {input.genes} \
			{params.flags}
		"""

### generate AllPredictions file
rule make_all_predictions:
	input: 
		predLists = expand(os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictions.tsv"), biosample=BIOSAMPLES_CONFIG['biosample'])
	params:
		output_dir = RESULTS_DIR
	conda:
		"../envs/abcenv.yml"
	output:
		allPred = os.path.join(RESULTS_DIR,"AllPredictions.txt.gz")
	resources:
		mem_gb=32,
		runtime_hr=6
	shell:
		"""			
		set +o pipefail;
		## make all predictions file 
		printf "chr\tstart\tend\tname\tTargetGene\tTargetGeneTSS\tCellType\tABC.Score\n" > {params.output_dir}/AllPredictions.txt
		for sample in {input.predLists}
		do
			cat $sample | sed 1d >> {params.output_dir}/AllPredictions.txt
		done
		gzip {params.output_dir}/AllPredictions.txt
		"""
