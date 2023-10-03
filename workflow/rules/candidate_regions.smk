rule make_candidate_regions:
	input:
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted"),
		accessibility = get_accessibility_file,
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", config['ref']['chrom_sizes'] + '.bed'),
	params:
		TSS = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'TSS'],
		chrom_sizes = config['ref']['chrom_sizes'],
		regions_blocklist = config['ref']['regions_blocklist'],
		peakExtendFromSummit = config['params_candidate']['peakExtendFromSummit'],
		nStrongestPeak = config['params_candidate']['nStrongestPeaks'],
		output_dir = os.path.join(RESULTS_DIR, "{biosample}", "Peaks"),
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/abcenv.yml"
	output: 
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
	resources:
		mem_mb=determine_mem_mb
	shell: 
		"""
		python {params.scripts_dir}/makeCandidateRegions.py \
			--narrowPeak {input.narrowPeak}\
			--bam {input.accessibility} \
			--outDir {params.output_dir} \
			--chrom_sizes {params.chrom_sizes} \
			--chrom_sizes_bed {input.chrom_sizes_bed} \
			--regions_blocklist {params.regions_blocklist} \
			--regions_includelist {params.TSS} \
			--peakExtendFromSummit {params.peakExtendFromSummit} \
			--nStrongestPeak {params.nStrongestPeak}
		"""
