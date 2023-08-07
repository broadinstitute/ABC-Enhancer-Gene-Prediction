## call macs2 -- if multiple accessibility inputs for one biosample, will aggregate into one output
rule call_macs_peaks: 
	input:
		accessibility = get_accessibility_file
	params:
		pval = config['params_macs']['pval'],
		genome_size = config['params_macs']['genome_size'],
		out_dir = config["predictions_results_dir"]
	conda:
		"../envs/abcenv.yml"
	output: 
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak")
	resources:
		mem_gb=32,
		runtime_hr=6
	shell: 
		""" 
		macs2 callpeak \
		-f AUTO \
		-g {params.genome_size} \
		-p {params.pval} \
		-n macs2 \
		--shift -75 \
		--extsize 150 \
		--nomodel \
		--keep-dup all \
		--call-summits \
		--outdir {params.out_dir}/{wildcards.biosample}/Peaks \
		-t {input.accessibility} 
		"""

rule gen_chrom_sizes_bed:
	input:
		chrom_sizes = config['chrom_sizes']
	output:
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", "chr_sizes.bed")
	resources:
		mem_gb=4,
		runtime_hr=1
	shell:
		"""
		awk 'BEGIN {{OFS="\t"}} {{if (NF > 0) print $1,"0",$2 ; else print $0}}' {input.chrom_sizes} > {output.chrom_sizes_bed}
		"""

## sort narrowPeaks
rule sort_narrowpeaks:
	input:
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak"),
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", "chr_sizes.bed")
	params:
		chrom_sizes = config['chrom_sizes']
	conda:
		"../envs/abcenv.yml"
	resources:
		mem_gb=4,
		runtime_hr=1
	output:
		narrowPeakSorted = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted")
	shell:
		"""
		# intersect to remove alternate chromosomes
		bedtools intersect -u -a {input.narrowPeak} -b {input.chrom_sizes_bed} | \
		bedtools sort -faidx {params.chrom_sizes} -i stdin > {output.narrowPeakSorted}
		"""
