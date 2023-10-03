class InvalidConfig(Exception):
	pass 

wildcard_constraints:
    threshold="\d+\.\d+",
	separator=".{0}|_",
	other_flags=".{0}|[^0-9]+"  # match empty strings or more flags

FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE = "threshold{threshold}{separator}{other_flags}"

MAX_MEM_MB = 250 * 1000  # 250GB

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	for _, file in input.items():
		if file.endswith('gz'):
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

def configure_hic_hashes(biosamples_config):
	"""
	Create a Map[hash(hic_dir), HIC_COLUMNS]
	This is done so that we don't recompute hic powerlaw fit
	when multiple biosamples have the same hic_info
	"""
	hic_hashes = {}
	hic_pairs = biosamples_config[HIC_COLUMNS].drop_duplicates()	
	for row in hic_pairs.values:
		hic_dir = row[0]
		if hic_dir:
			# Map only contains values with hic_directories
			hic_hashes[get_hic_dir_hash(row)] = row
	return hic_hashes

def get_accessibility_file(wildcards):
	# Inputs have been validated so only DHS or ATAC is provided
	return BIOSAMPLES_CONFIG.loc[wildcards.biosample, "DHS"] or BIOSAMPLES_CONFIG.loc[wildcards.biosample, "ATAC"]

def get_hic_dir_hash(hic_info_row):
	return hashlib.sha1(str(hic_info_row).encode()).hexdigest()[:8]


def get_hic_powerlaw_fit_file(wildcards):
	return os.path.join(_get_hic_powerlaw_fit_dir(wildcards), "hic.powerlaw.tsv")

def _get_hic_powerlaw_fit_dir(wildcards):
	"""
	If HiC is provided, we store the fit in the HiC hash folder. Otherwise
	we store under the biosamples folder
	"""
	row = BIOSAMPLES_CONFIG.loc[wildcards.biosample, HIC_COLUMNS].values
	hic_dir = row[0]
	if hic_dir:
		return os.path.join(RESULTS_DIR, "HiC_Powerlaw", get_hic_dir_hash(row))
	else:
		return os.path.join(RESULTS_DIR, "HiC_Powerlaw", wildcards.biosample)

def _validate_accessibility_feature(row: pd.Series):
	if row["DHS"] and row["ATAC"]:
		raise InvalidConfig("Can only specify one of DHS or ATAC for accessibility")
	if not (row["DHS"] or row["ATAC"]):
		raise InvalidConfig("Must provide either DHS or ATAC accessibility file")

def _validate_hic_info(row: pd.Series):
	if row["HiC_dir"]:
		if not (row["HiC_type"] and row["HiC_resolution"]):
			raise InvalidConfig("Must provide HiC type and resolution with directory")
	else:
		if not (row["HiC_gamma"] and row["HiC_scale"]):
			raise InvalidConfig(
				"If HiC dir not provided, you must provide HiC gamma and scale "
				"for powerlaw estimate for contact"
			)

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