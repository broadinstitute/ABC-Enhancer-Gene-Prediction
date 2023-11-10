
rule generate_qc_plot_and_summary:
	input: 
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
		neighborhoodDirectory = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods"),
		enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		chrom_sizes = config['ref']['chrom_sizes'],
	params:
		output_dir = os.path.join(RESULTS_DIR, "{biosample}", "Metrics"),
		scripts_dir = SCRIPTS_DIR,
		gamma = config['params_predict']['hic_gamma'],
		scale = config['params_predict']['hic_scale'],
	conda:
		"../envs/abcenv.yml"
	output:
		qc_summary = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", f"QCSummary_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		qc_plots = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", f"QCPlots_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.pdf")
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
		python {params.scripts_dir}/grabMetrics.py \
			--outdir {params.output_dir} \
			--output_qc_summary {output.qc_summary} \
			--output_qc_plots {output.qc_plots} \
			--macs_peaks {input.candidateRegions} \
			--neighborhood_outdir {input.neighborhoodDirectory} \
			--preds_file {input.enhPredictionsFull} \
			--chrom_sizes {input.chrom_sizes} \
			--hic_gamma {params.gamma} \
			--hic_scale {params.scale} 
		"""
