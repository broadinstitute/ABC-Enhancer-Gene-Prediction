
rule generate_qc_plot_and_summary:
	input: 
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
		neighborhoodDirectory = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods"),
		enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", f"EnhancerPredictionsFull_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		chrom_sizes = config['chrom_sizes'],
		powerlaw_params_tsv = get_hic_powerlaw_fit_file
	params:
		output_dir = os.path.join(RESULTS_DIR, "{biosample}", "Metrics")
	conda:
		"../envs/abcenv.yml"
	output:
		qc_summary = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", f"QCSummary_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.tsv"),
		qc_plots = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", f"QCPlots_{FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE}.pdf")
	resources:
		mem_gb=32,
		runtime_hr=6
	shell:
		"""
		python workflow/scripts/grabMetrics.py \
			--outdir {params.output_dir} \
			--output_qc_summary {output.qc_summary} \
			--output_qc_plots {output.qc_plots} \
			--macs_peaks {input.candidateRegions} \
			--neighborhood_outdir {input.neighborhoodDirectory} \
			--preds_file {input.enhPredictionsFull} \
			--chrom_sizes {input.chrom_sizes} \
			--powerlaw_params_tsv {input.powerlaw_params_tsv}
		"""
