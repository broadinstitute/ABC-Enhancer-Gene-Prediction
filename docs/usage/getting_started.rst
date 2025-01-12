.. _GettingStarted:

Getting Started
===============

Running ABC only requires an accessibility file, either DNase-seq (bam) or ATAC-seq (tagAlign)

Installation
------------

#. Download the repo locally from github
	- ``git clone git@github.com:broadinstitute/ABC-Enhancer-Gene-Prediction.git``
#. Make sure you have conda & mamba installed
	- `<https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
	- Make sure you're not using strict channel priorities: ``conda config --set channel_priority flexible``. Otherwise, you may encounter package conflicts later when installing abc. 
	- To install mamba: ``conda create -n mamba -c conda-forge mamba -y``
		- We recommend mamba as using conda can take 1hr+ for setup
		- See troubleshooting page if you run into issues


Setup Conda Environment
-----------------------
Creating the abc conda environment may take a while (~15min with mamba. > 1hr with conda)

.. code-block:: console

	$ conda activate mamba
	$ mamba env create -f workflow/envs/abcenv.yml
	$ conda activate abc-env


- Use `workflow/envs/release.yml` if you'd like the exact same abc conda environment as when ABC was released (we recommend using the abcenv.yml for most cases)

Running ABC
-----------
By default, the snakemake config file will run a test of ABC using data from K562 cells on chromosome 22 using reference genome hg38.

.. code-block:: console

	$ (abc-env) atan5133@NGPFWJVWWQ ABC-Enhancer-Gene-Prediction % snakemake -j1

You can see what commands snakemake will run

.. code-block:: console

	$ (abc-env) atan5133@NGPFWJVWWQ ABC-Enhancer-Gene-Prediction % snakemake -n -p

To read about each step, check out :ref:`ABC-methods`

The predictions will be stored in the ``{ABC_DIR}/results/{biosample_name}/Predictions`` folder. 
``EnhancerPredictions_threshold_.*.tsv`` contains the predicted ABC E-G links that meet the ABC Score threshold.
``EnhancerPredictionsAllPutative.tsv.gz`` contains all (unthresholded) E-G links with the ABC Score.

To sanity check your output from ABC, you can check out the QC metrics in the ``{ABC_DIR}/results/{biosample_name}/Metrics`` folder. 
For comparison, you can find the QC plots for our K562 run `here <https://drive.google.com/file/d/1fyd7ONKDgP646fOIafJhXcXnAk_6LCi1/view?usp=sharing>`_.
The metrics includes plots of things such as number of enhancers per gene and number of enhancer-genes per chromosome.


Configuring ABC
---------------

The primary configuration file for ABC is `config/config.yaml
<https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/config/config.yaml>`_


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

The rule specific params are explained in the :ref:`ABC-methods` section.

Genome Builds
-------------
The default reference file params in the config.yaml file are programmed for hg38 genome. To use a different genome, change the reference files and specify the genomize size parameter under `params_macs`.


BiosampleTable Specifications
-----------------------------
`chr22 example <https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/config/config_biosamples_chr22.tsv>`_

biosamples config is a tsv separated file with the following columns

#. Biosample 
	- Name to associate with your sample. e.g K562
#. DHS
	- DNAse-seq BAM file (sorted w/ .bai index file existence)
	- Can pass in multiple files separated by ``,``
#. ATAC
	- Bulk or single cell ATAC-seq TagAlign file (sorted w/ Tabix .tbi index file existence)
	- Can pass in multiple files separated by ``,``
#. H3K27ac
	- H3K27ac ChIP seq BAM file (sorted w/ .bai index file existence)
	- Can pass in multiple files separated by ``,``
#. default_accessibility_feature
	- Choices: "DHS", "ATAC" (If you provided DHS BAM file, you would put "DHS" here)
#. HiC_file
	- Filepath/link to a .hic file (recommended) or hic directory for the biosample cell type. 
	- If not provided, powerlaw is used to approximate contact
	- Examples: 
		- if filepath: `/path/to/k562.hic`
		- if link: `https://www.encodeproject.org/files/ENCFF621AIY/@@download/ENCFF621AIY.hic`
		- if directory: `/path/to/HiC`
#. HiC_type
	- Choices: hic, juicebox, avg, bedpe
	- If you passed in a .hic file, use ``hic``
	- If you dumped hic into a directory via JuicerTools, use ``juicebox``
	- If you have a bedpe file for contact, it should be a tab delimited file containing 8 columns (chr1,start1,end1,chr2,start2,end2,name,score)
#. HiC_resolution (int)
	- Currently only 5KB (kilobases) is supported
	- 5KB means dna regions are bucketed into 5KB bins and we measure contact between those bins
#. alt_TSS (optional; not recommended to fill)
	- Alternative TSS reference file 
#. alt_genes (optional; not recommended to fill)
	- Alternative Gene bound reference file

Required columns
	- biosample
	- DHS or ATAC file
	- default_accessibility_feature

If you don't have any cell specific HiC data, the recommendation is to not fill in any of the HiC columns, which will 
lead to using the powerlaw as the contact metric.


There is validation in Snakemake to make sure you provide the required inputs when running. 
The rest of the columns are optional, but providing them may help improve prediction performance.

You can run ABC on multiple biosamples via multiple rows in the tsv file. 


