rule create_neighborhoods:
	input:		
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", "chr_sizes.bed")
	params:
		DHS = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'DHS'] if BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'DHS'] else '',
		ATAC = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'ATAC'] if BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'ATAC'] else '',
		default = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'default_accessibility_feature'],
		H3K27ac = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'H3K27ac'] if BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'H3K27ac'] else '',
		genes = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'genes'],
		ubiquitous_genes = config['ubiquitous_genes'],
		chrom_sizes = config['chrom_sizes'],
		qnorm = config['params_neighborhoods']['qnorm'],
	conda:
		"../envs/abcenv.yml"
	output: 
		enhList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
		geneList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
		neighborhoodDirectory = directory(os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods")),
	shell:
		"""
		# get sorted & unique gene list
		# intersect first to remove alternate chromosomes
		bedtools intersect -u -a {params.genes} -b {input.chrom_sizes_bed} | \
		bedtools sort -faidx {params.chrom_sizes} -i stdin | \
		uniq > {params.genes}.sorted.uniq
						
		python workflow/scripts/run.neighborhoods.py \
			--candidate_enhancer_regions {input.candidateRegions} \
			--DHS {params.DHS} \
			--ATAC {params.ATAC} \
			--default_accessibility_feature {params.default} \
			--chrom_sizes {params.chrom_sizes} \
			--outdir {output.neighborhoodDirectory} \
			--genes {params.genes}.sorted.uniq \
			--ubiquitously_expressed_genes {params.ubiquitous_genes} \
			--qnorm {params.qnorm} \
			--H3K27ac {params.H3K27ac}
		"""