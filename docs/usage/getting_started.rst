Getting Started
===============

Installation
------------

#. Download the repo locally from github
	- ``git clone git@github.com:broadinstitute/ABC-Enhancer-Gene-Prediction.git``
	- Utilize the **dev** branch: ``git checkout dev``
#. Make sure you have conda & mamba installed
	- `<https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
	- To install mamba: ``conda install -n base -c conda-forge mamba``
		- We recommend mamba as using conda can take 1hr+ for setup


Setup Conda Environment
-----------------------
Creating the abc conda environment may take a while (~15min with mamba. > 1hr with conda)

.. code-block:: console

	$ mamba env create -f workflow/envs/abcenv.yml
	$ conda activate abc-env

- Use `workflow/envs/release.yml` if you'd like the same exact conda environment

Running ABC
-----------
By default, the snakemake config file will run a test of ABC using data from K562 cells on chromosome 22 using reference genome hg38.

.. code-block:: console

	$ (abc-env) atan5133@NGPFWJVWWQ ABC-Enhancer-Gene-Prediction % snakemake -j1

You can see what commands snakemake will run

.. code-block:: console

	$ (abc-env) atan5133@NGPFWJVWWQ ABC-Enhancer-Gene-Prediction % snakemake -n -p

To read about each step, check out :ref:`ABC-methods`


Configuring ABC
---------------

The primary configuration file for ABC is `config/config.yaml
<https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/dev/config/config.yaml>`_


*First couple lines in config/config.yaml*

.. code-block::

	### INPUT DATA
	# Add your inputs here
	biosamplesTable: "config/config-biosamples-chr22.tsv" 

To run ABC with your own specified data, create a **config-biosamples.tsv** file and replace ``"config/config-biosamples-chr22.tsv"`` with your file name. See section below for more info on how to create the biosample tsv file

Reference files
	- chrom_sizes: chromosome sizes file
		- FORMAT: TSV with 2 columns: chromsome (str), size (int) 
	- regions_blocklist: enhancer/promoter sequences to exclude from the model
		- FORMAT: BED 
	- ubiquitous_genes: genes that are always expressed, regardless of cell type (these genes do not typically have distal enhancers and so are flagged by the pipeline)
		- FORMAT: TSV with 1 column: gene name (str)
	- genes: list of genes corresponding with the genome
		- FORMAT: BED6 with ENSEMBL_ID as 7th column 
	- genome_tss: 500bp TSS region for each gene in the genes file
		- FORMAT: BED6 with ENSEMBL_ID as 7th column 



The default reference file params in the config.yaml file are programmed for hg38 genome. To use a different genome, change the reference files and specify the genomize size parameter under `params_macs`.

The rule specific params are explained in the :ref:`ABC-methods` section.



BiosampleTable Specifications
-----------------------------
`chr22 example <https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/dev/config/config_biosamples_chr22.tsv>`_

biosamples config is a tsv separated file with the following columns

#. Biosample 
	- Name to associate with your sample. e.g K562
#. DHS
	- DNAse-seq BAM file (sorted w/ .bai index file existence)
#. ATAC
	- ATAC-seq TagAlign file (sorted w/ Tabix .tbi index file existence)
#. H3K27ac
	- H3K27ac ChIP seq BAM file (sorted w/ .bai index file existence)
#. default_accessibility_feature
	- Choice: "DHS", "ATAC" (If you provided DHS BAM file, you would put "DHS" here)
#. HiC_dir
	- HiC directory for the biosample cell type. If not provided, powerlaw is used to approximate contact
#. HiC_type
	- e.g juicebox, avg, bed   (*explain what this means*)
#. HiC_resolution (int)
	- resolution of the HiC data  (*explain what this means*)
#. HiC_gamma (float)
	- represents the HiC data powerlaw fit slope   (*explain what this means*)
#. HiC_scale (float)
	- represents the HiC data powerlaw fit intercept   (*explain what this means*)
#. alt_TSS (optional; not recommended to fill)
	- Alternative TSS reference file 
#. alt_genes (optional; not recommended to fill)
	- Alternative Gene bound reference file

Required columns
	- biosample
	- DHS or ATAC
	- default_accessibility_feature
	- HiC info (dir, type, resolution) or powerlaw params (HiC gamma and scale)

There is validation in Snakemake to make sure you provide the required inputs when running. 
The rest of the columns are optional, but providing them may help improve prediction performance.

You can run ABC on multiple biosamples via multiple rows in the tsv file. 


