def _get_run_predictions_hic_params(wildcards):
	hic_dir = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_dir"]
	hic_type = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_type"]
	hic_resolution = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_resolution"]
	if hic_dir:
		return f"--hic_dir {hic_dir} --hic_type {hic_type} --hic_resolution {hic_resolution}"
	else:
		return "--score_column powerlaw.Score"

rule create_predictions:
	input:
		enhancers = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
		genes = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
		powerlaw_params_tsv = get_hic_powerlaw_fit_file
	params:
		cellType = lambda wildcards: wildcards.biosample, 
		output_dir = lambda wildcards: os.path.join(RESULTS_DIR, wildcards.biosample, "Predictions"),
		score_column = config['params_filter_predictions']['score_column'],
		hic_params = _get_run_predictions_hic_params,
		chrom_sizes = config['chrom_sizes'],
		flags = config['params_predict']['flags'],
		accessibility_feature = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'default_accessibility_feature'],
	conda:
		"../envs/abcenv.yml"
	output: 
		allPutative = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		allPutativeNonExpressed = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"),
	resources:
		mem_gb=128,
		runtime_hr=6
	shell:
		"""
		python workflow/scripts/predict.py \
			--enhancers {input.enhancers} \
			--outdir {params.output_dir} \
			--powerlaw_params_tsv {input.powerlaw_params_tsv} \
			--score_column {params.score_column} \
			--chrom_sizes {params.chrom_sizes} \
			--accessibility_feature {params.accessibility_feature} \
			--cellType {params.cellType} \
			--genes {input.genes} \
			{params.hic_params} \
			{params.flags}
		"""

rule filter_predictions:
	input: 
		allPutative = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		allPutativeNonExpressed = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"),
	params:
		score_column = config['params_filter_predictions']['score_column'],
		threshold = config['params_filter_predictions']['threshold'],
		include_self_promoter = config['params_filter_predictions']['include_self_promoter'],
		only_expressed_genes = config['params_filter_predictions']['only_expressed_genes']
	conda:
		"../envs/abcenv.yml"
	output:
		enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		enhPredictionsFullBed = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.bed"),
		genePredictionsStats = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"GenePredictionStats_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
	resources:
		mem_gb=32,
		runtime_hr=6
	shell:
		"""
		python workflow/scripts/filter_predictions.py \
			--output_tsv_file {output.enhPredictionsFull} \
			--output_bed_file {output.enhPredictionsFullBed} \
			--output_gene_stats_file {output.genePredictionsStats} \
			--pred_file {input.allPutative} \
			--pred_nonexpressed_file {input.allPutativeNonExpressed} \
			--score_column {params.score_column} \
			--threshold {params.threshold} \
			--include_self_promoter {params.include_self_promoter} \
			--only_expressed_genes {params.only_expressed_genes}
		"""
