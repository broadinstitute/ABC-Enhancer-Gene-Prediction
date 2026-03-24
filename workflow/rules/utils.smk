import pandas as pd
import time
import random
from pandas.errors import EmptyDataError

class InvalidConfig(Exception):
	pass 

wildcard_constraints:
	threshold=r"\d+\.\d+",
	separator=r".{0}|_",
	other_flags=r".{0}|[^0-9]+"  # match empty strings or more flags

FILTERED_PREDICTION_FILE_FORMAT_TEMPLATE = "threshold{threshold}{separator}{other_flags}"
DEFAULT_THRESHOLD = .02

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

def make_paths_absolute(obj, base_path):
	"""
	Use absolute paths to be compatible with github submodules
	Recursively go through the dictionary and convert relative paths to absolute paths.
	"""
	if isinstance(obj, dict):
		for key, value in obj.items():
			obj[key] = make_paths_absolute(value, base_path)
	elif isinstance(obj, str):
		# We assume all strings are paths. If converting the string
		# to an absolute path results in a valid file, then the str was a path
		new_file = os.path.join(base_path, obj)
		if os.path.exists(new_file):
			return new_file
	return obj

def determine_threshold(biosample):
	# config takes priority
	config_threshold = config["params_filter_predictions"]["threshold"]
	if config_threshold:
		return config_threshold
	biosample_row = BIOSAMPLES_CONFIG[BIOSAMPLES_CONFIG["biosample"] == biosample].iloc[0]
	hic_type = biosample_row["HiC_type"]
	if hic_type == None:
		hic_type = "powerlaw"
	elif hic_type == "avg":
		hic_type = "avg_hic"
	elif hic_type == "hic":
		hic_type = "intact_hic"
	matching_row = ABC_THRESHOLDS[
        (ABC_THRESHOLDS["accessibility"] == biosample_row["default_accessibility_feature"])
        & (ABC_THRESHOLDS["has_h3k27ac"] == bool(biosample_row["H3K27ac"]))
        & (ABC_THRESHOLDS["hic_type"] == hic_type)
    ]
	if len(matching_row) == 0:
		print(f"Threshold not found for biosample: {biosample}. Using default threshold of {DEFAULT_THRESHOLD}")
		threshold = DEFAULT_THRESHOLD
	else:
		threshold = matching_row.iloc[0]["threshold"]
	return threshold

def determine_filtered_prediction_file_format(threshold, config):
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

def enable_retry(func, func_args={}, max_attempts=3, delay=0.5):
	"""
	To prevent EmptyDataError race condition when using SLURM ro launch jobs as processes
	Assuming the EmptyDataError is caused by a file caching or synchronization lag
	Retry with delay

	@Param
	func:  Function to retry
	func_args:  Dictionary of kwargs for function
	max_attempts:  Maximum number of attempts allowable before raising error
	delay: minimum delay before retry
	"""
	for attempt in range(max_attempts):
		try:
			return func(**func_args)
		except Exception as e:
			if attempt == max_attempts - 1:
				raise
			sleep_time = delay + random.uniform(0, 0.5)
			time.sleep(sleep_time)
	return None

def load_biosamples_config(config, validate_inputs_exist=True):
	biosamples_config = enable_retry(
		pd.read_csv, 
		func_args={'filepath_or_buffer': config["biosamplesTable"], 'sep': "\t"}
	).replace([np.nan], [None]).infer_objects(copy=False).set_index("biosample", drop=False)
	biosamples_config["HiC_resolution"] = biosamples_config["HiC_resolution"].fillna(0).astype(int)
	_validate_biosamples_config(biosamples_config, validate_inputs_exist)
	_configure_tss_and_gene_files(biosamples_config)
	return biosamples_config

def load_abc_thresholds(config):
	file = config["ref"]["abc_thresholds"]
	return pd.read_csv(file, sep='\t')

def get_accessibility_files(wildcards):
	# Inputs have been validated so only DHS or ATAC is provided
	files = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "DHS"] or BIOSAMPLES_CONFIG.loc[wildcards.biosample, "ATAC"]
	return files.split(",")

def get_activity_files(wildcards):
	# for neighborhoods step, to trigger download of necessary inputs
	files = get_accessibility_files(wildcards)
	k27ac_value = BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K27ac"]
	if k27ac_value:
		k27ac_files = k27ac_value.split(",")
		files = files + k27ac_files
	return files

def get_index_file(filepath):
	"""
	Return the expected index file path for a given input file.
	- BAM files (.bam) -> .bam.bai
	- tagAlign files (.gz) -> .gz.tbi
	- BigWig files -> None (no index needed)
	"""
	if filepath.endswith(".bam"):
		return filepath + ".bai"
	elif "tagAlign" in filepath and filepath.endswith(".gz"):
		return filepath + ".tbi"
	# BigWig and other files don't need indexes
	return None

def get_index_files_for_list(filepaths):
	"""Return list of index files for files that need them."""
	indexes = []
	for f in filepaths:
		idx = get_index_file(f)
		if idx:
			indexes.append(idx)
	return indexes

def get_accessibility_index_files(wildcards):
	"""Return index files for accessibility BAM/tagAlign files."""
	files = get_accessibility_files(wildcards)
	return get_index_files_for_list(files)

def get_activity_index_files(wildcards):
	"""Return index files for all activity files (accessibility + H3K27ac)."""
	files = get_activity_files(wildcards)
	return get_index_files_for_list(files)

def _validate_accessibility_feature(row: pd.Series):
	if row["DHS"] and row["ATAC"]:
		raise InvalidConfig("Can only specify one of DHS or ATAC for accessibility")
	if not (row["DHS"] or row["ATAC"]):
		raise InvalidConfig("Must provide either DHS or ATAC accessibility file")

def _validate_hic_info(row: pd.Series):
	if row["HiC_file"]:
		if not (row["HiC_type"] and row["HiC_resolution"]):
			raise InvalidConfig("Must provide HiC type and resolution with file")
		if row["HiC_resolution"] != 5000:
			raise InvalidConfig("Only 5kb resolution supported at the moment")

def _is_url(path):
	"""Check if a path is a URL."""
	return path.startswith("http://") or path.startswith("https://")

def _validate_input_files_exist(row: pd.Series):
	"""
	Validate that input files and their indexes exist.
	Raises InvalidConfig with informative message if files are missing.
	Skips validation for URLs (e.g., HiC files from ENCODE).
	"""
	biosample = row["biosample"]
	missing_files = []
	missing_indexes = []

	# Check accessibility files (DHS or ATAC)
	access_files_str = row["DHS"] or row["ATAC"]
	if access_files_str:
		for f in access_files_str.split(","):
			f = f.strip()
			if not os.path.exists(f):
				missing_files.append(f)
			else:
				idx = get_index_file(f)
				if idx and not os.path.exists(idx):
					missing_indexes.append(idx)

	# Check H3K27ac files
	if row["H3K27ac"]:
		for f in row["H3K27ac"].split(","):
			f = f.strip()
			if not os.path.exists(f):
				missing_files.append(f)
			else:
				idx = get_index_file(f)
				if idx and not os.path.exists(idx):
					missing_indexes.append(idx)

	# Check HiC file (skip URLs - they are fetched at runtime)
	if row["HiC_file"] and not _is_url(row["HiC_file"]) and not os.path.exists(row["HiC_file"]):
		missing_files.append(row["HiC_file"])

	# Report errors
	if missing_files or missing_indexes:
		error_msg = f"Input file validation failed for biosample '{biosample}':\n"
		if missing_files:
			error_msg += "\n  Missing files:\n"
			for f in missing_files:
				error_msg += f"    - {f}\n"
		if missing_indexes:
			error_msg += "\n  Missing index files:\n"
			for f in missing_indexes:
				error_msg += f"    - {f}\n"
			error_msg += "\n  To create missing indexes:\n"
			error_msg += "    - For BAM files: samtools index <file.bam>\n"
			error_msg += "    - For tagAlign.gz: tabix -p bed <file.tagAlign.gz>\n"
		raise InvalidConfig(error_msg)

def _validate_biosamples_config(biosamples_config, validate_inputs_exist=True):
	"""
	Throw exception if a row needs to be fixed
	"""
	for _, row in biosamples_config.iterrows():
		_validate_hic_info(row)
		_validate_accessibility_feature(row)
		if validate_inputs_exist:
			_validate_input_files_exist(row)

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
