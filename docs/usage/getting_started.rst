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
By default, ABC is configured to run on a a K562 cell type for chromosome 22 using reference genome hg38.

.. code-block:: console

	# (abc-env) atan5133@NGPFWJVWWQ ABC-Enhancer-Gene-Prediction % snakemake -j1

To see the multiple steps that snakemake performs, check out :ref:`ABC-methods`

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

The other params in the config.yaml file are programmed for hg38 with corresponding reference files. To use a different genome, change the reference file and specify the genomize size parameter under `params_macs`.

BiosampleTable Specifications
-----------------------------
`chr22 example <https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/dev/config/config-biosamples-chr22.tsv>`_

biosamples config is a tsv separated file with the following columns

#. Biosample 
	- Name to associate with your sample. e.g K562
#. DHS
	- DNAse-seq BAM file (indexed and sorted)
#. ATAC
	- ATAC-seq BAM file (indexed and sorted)
#. H3K27ac
	- H3K27ac ChIP seq BAM file
#. default_accessibility_feature
	- Either DHS or ATAC
#. HiC_dir
	- HiC directory for the biosample cell type. If not provided, uses avg hi-c
#. HiC_type
	- e.g juicebox, avg, bed
#. HiC_resolution (int)
	- resolution of the HiC data
#. HiC_gamma (float)
	- represents the HiC data powerlaw fit slope
#. HiC_scale (float)
	- represents the HiC data powerlaw fit intercept
#. alt_TSS
	- Alternative TSS reference file 
#. alt_genes
	- Alternative Gene bound reference file

Required columns
- biosample
- DHS or ATAC
- default_accessibility_feature
- either HiC info (dir, type, resolution) and/or HiC gamma and scale

There is validation in Snakemake to make sure you provide the required inputs when running. 

The rest of the columns are optional, but providing them may help improve prediction performance.

You can run ABC on multiple biosamples by inputting multiple rows in the tsv file. 


