def _get_write_powerlaw_params(wildcards) -> str:
	if wildcards.KEY in HIC_HASHES:
		hic_dir, hic_type, hic_resolution, hic_gamma, hic_scale = HIC_HASHES[wildcards.KEY]
		params = f"--hic_dir {hic_dir} --hic_type {hic_type} --hic_resolution {hic_resolution}"
		if hic_gamma and hic_scale:
			params += f" --hic_gamma {hic_gamma} --hic_scale {hic_scale}"
		return params
	else:
		# Otherwise, KEY is a biosample and there's no HiC directory associated
		hic_gamma = BIOSAMPLES_CONFIG.loc[wildcards.KEY, "HiC_gamma"]
		hic_scale = BIOSAMPLES_CONFIG.loc[wildcards.KEY, "HiC_scale"]
		# B/c of config file validation, we know gamma/scale values exist here
		return f"--hic_gamma {hic_gamma} --hic_scale {hic_scale}"
	
rule write_powerlaw_params:
	params:
		params = _get_write_powerlaw_params,
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/abcenv.yml"
	output:
		powerlaw_params_tsv = os.path.join(RESULTS_DIR, "HiC_Powerlaw", "{KEY}", "hic.powerlaw.tsv")
	resources:
		mem_mb=determine_mem_mb
	shell:
		"""
		python {params.scripts_dir}/write_powerlaw_params.py \
			{params.params} \
			--output_file {output.powerlaw_params_tsv}
		"""