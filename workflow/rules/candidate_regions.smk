rule make_candidate_regions:
	input:
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted"),
	params:
		accessibility = get_accessibility_file,
		TSS = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'TSS'],
		chrom_sizes = config['chrom_sizes'],
		regions_blocklist = config['regions_blocklist'],
		peakExtendFromSummit = config['params_candidate']['peakExtendFromSummit'],
		nStrongestPeak = config['params_candidate']['nStrongestPeaks'],
		output_dir = RESULTS_DIR
	conda:
		"../envs/abcenv.yml"
	output: 
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
	shell: 
		"""
		python workflow/scripts/makeCandidateRegions.py \
			--narrowPeak {input.narrowPeak}\
			--bam {params.accessibility} \
			--outDir {params.output_dir}/{wildcards.biosample}/Peaks \
			--chrom_sizes {params.chrom_sizes} \
			--regions_blocklist {params.regions_blocklist} \
			--regions_includelist {params.TSS} \
			--peakExtendFromSummit {params.peakExtendFromSummit} \
			--nStrongestPeak {params.nStrongestPeak}
		"""