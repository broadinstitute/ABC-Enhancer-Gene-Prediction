rule create_neighborhoods:
	input:		
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
	params:
		DHS = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "DHS"] or '',
		ATAC = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "ATAC"] or '',
		default = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'default_accessibility_feature'],
		H3K27ac = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K27ac"] or '',
		genes = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'genes'],
		ubiquitous_genes = config['ref']['ubiquitous_genes'],
		chrom_sizes = config['ref']['chrom_sizes'],
		qnorm = f"--qnorm {config['ref']['qnorm']}" if config['params_neighborhoods']['use_qnorm'] else "",
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/abcenv.yml"
	output: 
		enhList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
		geneList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
		neighborhoodDirectory = directory(os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods")),
		processed_genes_file = os.path.join(RESULTS_DIR, "{biosample}", "processed_genes_file.bed"),
	resources:
		mem_mb=32*1000
	shell:
		"""
		# get sorted & unique gene list
		# intersect first to remove alternate chromosomes
		bedtools intersect -u -a {params.genes} -b {input.chrom_sizes_bed} | \
		bedtools sort -faidx {params.chrom_sizes} -i stdin | \
		uniq > {output.processed_genes_file}
						
		python {params.scripts_dir}/run.neighborhoods.py \
			--candidate_enhancer_regions {input.candidateRegions} \
			--DHS {params.DHS} \
			--ATAC {params.ATAC} \
			--default_accessibility_feature {params.default} \
			--chrom_sizes {params.chrom_sizes} \
			--chrom_sizes_bed {input.chrom_sizes_bed} \
			--outdir {output.neighborhoodDirectory} \
			--genes {output.processed_genes_file} \
			--ubiquitously_expressed_genes {params.ubiquitous_genes} \
			--H3K27ac {params.H3K27ac} \
			{params.qnorm} 
		"""