## call macs2 -- if multiple accessibility inputs for one biosample, will aggregate into one output
rule call_macs_peaks: 
	input:
		accessibility = get_accessibility_file
	params:
		pval = config['params_macs']['pval'],
		out_dir = config["predictions_results_dir"]
	conda:
		"../envs/abcenv.yml"
	output: 
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak")
	shell: 
		""" 
		macs2 callpeak -f AUTO -g hs -p {params.pval} -n macs2 --call-summits --outdir {params.out_dir}/{wildcards.biosample}/Peaks -t {input.accessibility}
		"""

## sort narrowPeaks
rule sort_narrowpeaks:
	input:
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak")
	params:
		chrom_sizes = config['chrom_sizes']
	conda:
		"../envs/abcenv.yml"
	output:
		narrowPeakSorted = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted")
	shell:
		"""
		# intersect first to remove alternate chromosomes
		bedtools intersect -u -a {input.narrowPeak} -b {params.chrom_sizes}.bed | \
		bedtools sort -faidx {params.chrom_sizes} -i stdin > {output.narrowPeakSorted}
		"""