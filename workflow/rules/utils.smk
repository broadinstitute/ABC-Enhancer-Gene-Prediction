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

MAX_MEM_MB = 250 * 1000  # 250GB

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

def load_biosamples_config(config):
	biosamples_config = enable_retry(
		pd.read_csv, 
		func_args={'filepath_or_buffer': config["biosamplesTable"], 'sep': "\t"}
	).replace([np.nan], [None]).set_index("biosample", drop=False)
	biosamples_config["HiC_resolution"] = biosamples_config["HiC_resolution"].replace([None], [0]).astype(int)
	_validate_biosamples_config(biosamples_config)
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
