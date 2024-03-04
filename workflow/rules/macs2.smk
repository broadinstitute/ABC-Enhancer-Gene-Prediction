## call macs2 -- if multiple accessibility inputs for one biosample, will aggregate into one output
rule call_macs_peaks: 
	input:
		accessibility = get_accessibility_files,
	params:
		pval = config['params_macs']['pval'],
		genome_size = config['params_macs']['genome_size'],
	conda:
		"../envs/abcenv.yml"
	output: 
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak")
	resources:
		mem_mb=determine_mem_mb
	shell: 
		"""
		if [[ "{input.accessibility}" == *tagAlign* ]]; then
			FORMAT="BED"
		else
			FORMAT="AUTO"
		fi

		macs2 callpeak \
		-f $FORMAT \
		-g {params.genome_size} \
		-p {params.pval} \
		-n macs2 \
		--shift -75 \
		--extsize 150 \
		--nomodel \
		--keep-dup all \
		--call-summits \
		--outdir {RESULTS_DIR}/{wildcards.biosample}/Peaks \
		-t {input.accessibility} 
		"""

rule generate_chrom_sizes_bed_file:
	input:
		chrom_sizes = config['ref']['chrom_sizes']
	output:
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
		awk 'BEGIN {{OFS="\t"}} {{if (NF > 0) print $1,"0",$2 ; else print $0}}' {input.chrom_sizes} > {output.chrom_sizes_bed}
		"""

## sort narrowPeaks
rule sort_narrowpeaks:
	input:
		narrowPeak = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak"),
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
	params:
		chrom_sizes = config['ref']['chrom_sizes']
	conda:
		"../envs/abcenv.yml"
	output:
		narrowPeakSorted = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted")
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
		# intersect to remove alternate chromosomes
		bedtools intersect -u -a {input.narrowPeak} -b {input.chrom_sizes_bed} | \
		bedtools sort -faidx {params.chrom_sizes} -i stdin > {output.narrowPeakSorted}
		"""
