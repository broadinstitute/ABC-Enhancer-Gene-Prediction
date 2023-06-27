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
Creating the abc conda environment may take a while (~30 min)

.. code-block:: console

	$ mamba env create -f workflow/envs/abcenv.yml
	$ conda activate abc-env

Running ABC
-----------
By default, ABC is configured to run on a a K562 cell type for chromosome 22 using reference genome hg38.

.. code-block:: console

	# (abc-env) atan5133@NGPFWJVWWQ ABC-Enhancer-Gene-Prediction % snakemake -j1

To see the multiple steps that snakemake performs, check out :ref:`ABC-steps`

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

The other params in the config.yaml file are programmed for hg38. If you wish to run against a different reference genome, you would have to provide the appropriate reference genome files. 

BiosampleTable Specifications
-----------------------------
`chr22 example <https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/dev/config/config-biosamples-chr22.tsv>`_

biosamples config is a tsv separate file with the following columns

#. Biosample 
	- Name to associate with your sample. e.g K562
#. DHS
	- DNAse-seq BAM file (indexed and sorted)
#. ATAC
	- ATAC-seq BAM file (indexed and sorted)
#. H3K27ac (Optional)
	- H3K27ac ChIP seq BAM file
#. default_accessibility_feature
	- Either DHS or ATAC
#. HiC_dir (optional)
	- HiC directory for the biosample cell type. If not provided, uses avg hi-c
#. HiC_type (optional)
#. alt_TSS (optional)
	- Alternative TSS reference file 
#. alt_genes (optional)
	- Alternative Gene bound reference file

Required columns

- biosample
- DHS or ATAC
- default_accessibility_feature

The rest of the columns are optional, but providing them may help improve prediction performance.

You can run ABC on multiple biosamples by inputting multiple rows.


