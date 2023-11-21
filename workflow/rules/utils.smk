class InvalidConfig(Exception):
	pass 

wildcard_constraints:
	threshold=r"\d+\.\d+",
	separator=r".{0}|_",
	other_flags=r".{0}|[^0-9]+"  # match empty strings or more flags

FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE = "threshold{threshold}{separator}{other_flags}"

MAX_MEM_MB = 250 * 1000  # 250GB

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 5  # assume gz compressesed the file <= 5x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

def convert_reference_files(config):
	# prefixes ABC path if config var is provided
	abs_path_prefix = config.get("ABC_DIR_PATH")
	if abs_path_prefix:
		for name, ref_file in config["ref"].items():
			config["ref"][name] = os.path.join(abs_path_prefix, ref_file)
	return config

def determine_filtered_prediction_file_format(config):
	threshold = config['params_filter_predictions']['threshold']
	include_self_promoter = config['params_filter_predictions']['include_self_promoter']
	only_expressed_genes = config['params_filter_predictions']['only_expressed_genes']
	if include_self_promoter or only_expressed_genes:
		separator = '_'
		other_flags = []
		if include_self_promoter:
			other_flags.append('self_promoter')
		if only_expressed_genes:
			other_flags.append('only_expr_genes')
		other_flags = "__".join(other_flags)
	else:
		separator = ''
		other_flags = ''
	return FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE.format(threshold=threshold, separator=separator, other_flags=other_flags)

def load_biosamples_config(config):
	biosamples_config = pd.read_csv(
		config["biosamplesTable"], sep="\t", na_values=""
	).replace([np.NaN], [None]).set_index("biosample", drop=False)
	biosamples_config["HiC_resolution"] = biosamples_config["HiC_resolution"].replace([None], [0]).astype(int)
	_validate_biosamples_config(biosamples_config)
	_configure_tss_and_gene_files(biosamples_config)
	return biosamples_config

def get_accessibility_files(wildcards):
	# Inputs have been validated so only DHS or ATAC is provided
	files = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "DHS"] or BIOSAMPLES_CONFIG.loc[wildcards.biosample, "ATAC"]
	return files.split(",")


def _validate_accessibility_feature(row: pd.Series):
	if row["DHS"] and row["ATAC"]:
		raise InvalidConfig("Can only specify one of DHS or ATAC for accessibility")
	if not (row["DHS"] or row["ATAC"]):
		raise InvalidConfig("Must provide either DHS or ATAC accessibility file")

def _validate_hic_info(row: pd.Series):
	if row["HiC_file"]:
		if not (row["HiC_type"] and row["HiC_resolution"]):
			raise InvalidConfig("Must provide HiC type and resolution with file")

def _validate_biosamples_config(biosamples_config):
	"""
	Throw exception if a row needs to be fixed
	"""
	for _, row in biosamples_config.iterrows():
		_validate_hic_info(row)
		_validate_accessibility_feature(row)

def _configure_tss_and_gene_files(biosamples_config):
	## get TSS and genefile names for each biosample 
	TSS_files = []
	gene_files = []
	for sample in biosamples_config['biosample']:
		tss_file = config['ref']['genome_tss']
		gene_file = config['ref']['genes']
		if biosamples_config.loc[sample, "alt_TSS"]:
			tss_file = biosamples_config.loc[sample, 'alt_TSS']
		if biosamples_config.loc[sample, "alt_genes"]:
			gene_file = biosamples_config.loc[sample, 'alt_genes']
		TSS_files.append(tss_file)
		gene_files.append(gene_file)
					
	biosamples_config["TSS"] = TSS_files
	biosamples_config["genes"] = gene_files
