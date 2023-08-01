rule gen_qc_plots:
	input: 
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
		neighborhoodDirectory = directory(os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods")),
		enhPredictionsFull = os.path.join(RESULTS_DIR, "{biosample}", "Predictions", "EnhancerPredictionsFull.tsv"),
		powerlaw_params_tsv = get_hic_powerlaw_fit_file
	params:
		output_dir = os.path.join(RESULTS_DIR, "{biosample}", "Metrics")
	conda:
		"../envs/abcenv.yml"
	output:
		qc_summary = os.path.join(RESULTS_DIR, "{biosample}", "Metrics", "QCSummary.p")
	shell:
		"""
		python workflow/scripts/grabMetrics.py \
			--outdir {params.output_dir} \
			--macs_peaks {input.candidateRegions} \
			--neighborhood_outdir {input.neighborhoodDirectory} \
			--preds_file {input.enhPredictionsFull} \
			--powerlaw_params_tsv {input.powerlaw_params_tsv}
		"""
