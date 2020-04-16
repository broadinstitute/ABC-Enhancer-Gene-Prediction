rule download_bam:
        message: "Executing est.sh to download files"
	        output: expand("bam_{assay}.tsv", assay=ASSAYS)
		        log: "logs/download.log"
			        shell: "bash est.sh"

				rule obtain_paired_single_files:
				        output: expand("{pairtype}_bam_{assays}.tsv", assays=ASSAYS, pairtype=endtype)
					        message: "Executing log_kristy.sh to download paired-end / single-end file information"
						        log: "logs/pe_se_download.log"
							        shell: "bash log_kristy.sh {wildcards.inputs}"
