from functools import partial

def _get_run_predictions_hic_params(wildcards):
	hic_file = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_file"]
	hic_type = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_type"]
	hic_resolution = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "HiC_resolution"]
	if hic_file:
		return f"--hic_file {hic_file} --hic_type {hic_type} --hic_resolution {hic_resolution}"
	else:
		return "--score_column powerlaw.Score"

rule create_predictions:
	input:
		enhancers = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
		genes = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
	params:
		cellType = lambda wildcards: wildcards.biosample, 
		output_dir = lambda wildcards: os.path.join(RESULTS_DIR, wildcards.biosample, "Predictions"),
		score_column = config['params_filter_predictions']['score_column'],
		hic_params = _get_run_predictions_hic_params,
		chrom_sizes = config['ref']['chrom_sizes'],
		flags = config['params_predict']['flags'],
		gamma = config['params_predict']['hic_gamma'],
		scale = config['params_predict']['hic_scale'],
		hic_pseudocount_distance = config['params_predict']['hic_pseudocount_distance'],
		accessibility_feature = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'default_accessibility_feature'],
		scripts_dir = SCRIPTS_DIR,
	conda:
		"../envs/abcenv.yml"
	output: 
		allPutative = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		allPutativeNonExpressed = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"),
	resources:
		mem_mb=partial(determine_mem_mb, min_gb=20)  # Use 100GB if using average HiC
	shell:
		"""
		python {params.scripts_dir}/predict.py \
			--enhancers {input.enhancers} \
			--outdir {params.output_dir} \
			--score_column {params.score_column} \
			--chrom_sizes {params.chrom_sizes} \
			--accessibility_feature {params.accessibility_feature} \
			--cellType {params.cellType} \
			--genes {input.genes} \
			--hic_gamma {params.gamma} \
			--hic_scale {params.scale} \
			--hic_pseudocount_distance {params.hic_pseudocount_distance} \
			{params.hic_params} \
			{params.flags}
		"""

rule combine_accessible_regions:
	input:
		prediction_files=expand(os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"), biosample=BIOSAMPLES_CONFIG["biosample"].to_list())
	output:
		combined_prediction_file=os.path.join(RESULTS_DIR, "combined_Predictions", "combined_accessible_peaks.tsv.gz"),
	params:
		biosamples=BIOSAMPLES_CONFIG["biosample"].to_list(),
		result_dir=os.path.join(RESULTS_DIR)
	resources:
		mem_mb=8*1000
	threads: 8
	shell:
		"""
		i=0
        for sample in {params.biosamples}
        do 
            if [ $i -eq 0 ]
            then
            	head -1 {params.result_dir}/$sample/Neighborhoods/EnhancerList.txt | cut -f1,2,3 |tr \'\\n\' \'\\t\' > {params.result_dir}/combined_Predictions/combined_accessible_peaks.tsv
				echo "celltype" >> {params.result_dir}/combined_Predictions/combined_accessible_peaks.tsv
				sed 1d {params.result_dir}/$sample/Neighborhoods/EnhancerList.txt | cut -f1,2,3 |sed "s/$/\t$sample/" >> {params.result_dir}/combined_Predictions/combined_accessible_peaks.tsv
            else
                sed 1d {params.result_dir}/$sample/Neighborhoods/EnhancerList.txt |cut -f1,2,3| sed "s/$/\t$sample/" >> {params.result_dir}/combined_Predictions/combined_accessible_peaks.tsv
            fi
            ((i=i+1))
        done
		cat {params.result_dir}/combined_Predictions/combined_accessible_peaks.tsv |pigz -p{threads} > {output.combined_prediction_file}
		rm {params.result_dir}/combined_Predictions/combined_accessible_peaks.tsv
		"""

rule combine_all_putative_predictions:
	input:
		all_putative_files=expand(os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"), biosample=BIOSAMPLES_CONFIG["biosample"].to_list())
	output:
		combined_all_putative_file=os.path.join(RESULTS_DIR, "combined_Predictions", "combined_only_expressed_genes_EnhancerPredictionsAllPutative.tsv.gz"),
	params:
		biosamples=BIOSAMPLES_CONFIG["biosample"].to_list(),
		result_dir=os.path.join(RESULTS_DIR)
	resources:
		mem_mb=8*1000
	threads: 16
	shell:
		"""
		touch {output.combined_all_putative_file}
		LC_ALL=C
        	for sample in {params.biosamples}
        	do 
            		zcat {params.result_dir}/$sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz | sed 1d |cut -f12,28|sort --parallel={threads} |uniq | pigz -p{threads} >> {output.combined_all_putative_file}
        	done
		"""
		#cat {params.result_dir}/combined_Predictions/combined_only_expressed_genes_EnhancerPredictionsAllPutative.tsv|pigz -p{threads} > {output.combined_all_putative_file}
		#rm {params.result_dir}/combined_Predictions/combined_only_expressed_genes_EnhancerPredictionsAllPutative.tsv

rule filter_predictions:
	input: 
		allPutative = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutative.tsv.gz"),
		allPutativeNonExpressed = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"),
	params:
		score_column = config['params_filter_predictions']['score_column'],
		threshold = lambda wildcards: determine_threshold(wildcards.biosample),
		include_self_promoter = config['params_filter_predictions']['include_self_promoter'],
		only_expressed_genes = config['params_filter_predictions']['only_expressed_genes'],
	conda:
		"../envs/abcenv.yml"
	output:
		enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		enhPredictionsFullBedpe = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.bedpe.gz"),
		enhPredictionsSlim = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictions_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		genePredictionsStats = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"GenePredictionStats_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv")
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
		python {SCRIPTS_DIR}/filter_predictions.py \
			--output_tsv_file {output.enhPredictionsFull} \
			--output_slim_tsv_file {output.enhPredictionsSlim} \
			--output_bed_file {output.enhPredictionsFullBedpe} \
			--output_gene_stats_file {output.genePredictionsStats} \
			--pred_file {input.allPutative} \
			--pred_nonexpressed_file {input.allPutativeNonExpressed} \
			--score_column {params.score_column} \
			--threshold {params.threshold} \
			--include_self_promoter {params.include_self_promoter} \
			--only_expressed_genes {params.only_expressed_genes}
		"""
