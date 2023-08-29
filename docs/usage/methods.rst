.. _ABC-methods:

Detailed Methods
================

To-do:  This page should include detailed methods for each of the main operations for ABC

We may want to include benchmarking figures or other methodological details e.g. from previous papers

Key concepts:

- Defining candidate elements
- Estimating enhancer activity
- Estimating enhancer-promoter 3D contact
- Making predictions with different combinations of input datasets
- Interpreting the ABC score

1. Defining candidate elements
------------------------------

'Candidate elements' are the set of putative enhancer elements for which ABC Scores will be computed. These include gene promoters. This step is encoded in the makeCnadidateRegions.py script.

Main inputs
	- narrowPeak file from MACS2 
	- DNase-seq of ATAC-seq input file

Output
	- Candidate regions in BED file format

Description:
	#. Count DNase-seq reads in each peak and retain the top N peaks (defined by --nStrongestPeak) with the most read counts
	#. Resize each of these N peaks to be a fixed number of base pairs centered on the peak summit (defined by --peakExtendFromSummit)
	#. Remove any regions listed in the 'blocklist' and include any regions listed in the 'includelist'
	#. Merge any overlapping regions


**Example:**

.. code-block:: console

	$ python workflow/scripts/makeCandidateRegions.py \
		--narrowPeak results/K562_chr22/Peaks/macs2_peaks.narrowPeak.sorted \
		--bam example_chr/chr22/ENCFF860XAE.chr22.sorted.se.bam \
		--outDir results/K562_chr22/Peaks \
		--chrom_sizes reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv \
		--regions_blocklist reference/hg38/GRCh38_unified_blacklist.bed \
		--regions_includelist example_chr/chr22/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.hg38.TSS500bp.bed \
		--peakExtendFromSummit 250 \
		--nStrongestPeak 150000

The method of defining candidate elements includes the following steps:

- Peak-calling with MACS2
- Resizing and merging regions
- Selecting the 150,000 strongest peaks (by read count)
- Adding promoters

1.1. Calling peaks with MACS2 [Rosa please edit]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Rosa to add details about the MACS2 methods calling


1.2. Resizing and merging regions [Rosa + Maya please edit]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
e.g. explain the logic:  We resize to 500 bp to count reads in and around peaks (esp. e.g. H3K27ac signal is surrounding the peak); peak callers are sensitive; 500 bp also a reasonable window for interpreting variants, since they're typically within this distance of a peak

Maya do you ahve a figure showing GWAS performance as a function of window size?


1.3. Selecting the top N peaks [Rosa please edit]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
explain the logic of this step, e.g. because experiments of different sequencing depths can have large differences in numbers of peaks; this affects ABC score via the denominator; so we include top 150K

  
1.4. Defining and adding gene promoters [Andreas please edit]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, we force the inclusion of gene promoters in the set of candidate elements. This is done because include promoters in the calculation; but sometimes the promoters of genes do not pass the threshold for the top 150,000 genes, which has a large effect on ABC due to the promoter receiving a high "3D contact" value in the ABC computation.

Note that the exact method of defining the promoter region for a gene has a strong influence on the ABC computation.

Describe how it is important how promoters are selected, and how changing the promoter list can impact ABC scores
	- First, the exact promoter used affects the ABC score for the gene corresponding to that promoter, because of 3D contacts (which can differ depending on the location of the promoter) and whether that promoter is in fact the dominant element used (the promoter is included as a candidate "enhancer" for itself, and contributes to the denominator of the ABC score).
	- Second, the promoter list used affects the ABC scores for other nearby genes, because force inclusion of these regions leads to more/larger regions being used which affects the ABC denominator.
	- Third, the promoter list used can affect downstream benchmarking analyses. For example, benchmarks that filter to just 'distal elements' that are not promoter might filter out elements called as promoters that are actually enhancers (e.g. promoters of lncRNAs that act as enhancers).

In practice, we provide a gene promoter file that we have used for various purposes that selects a single canonical promoter per gene. 
	- describe provenance of the gene promoter file(s) including in ABC repo (for human and mouse)
	- Changing the promoter for a single gene, e.g. to accommodate a specific alternative transcription start site of a gene of interest, is likely to be okay and not globally affect predictions
	- However, caution is warranting in making more extensive changes to the promoter list. Note again that including a much larger promoter list, e.g. including lncRNAs or including all possible transcription start sites for all isoforms for a gene, is likely to change the global properties of the ABC score and is not recommended without calibration of scores (see section on Interpreting the ABC score below)




2. Estimating enhancer activity
-------------------------------

In this step, we estimate the 'enhancer activity' of candidate elements by counting reads from ATAC, DNase-seq, and/or H3K27ac ChIP-seq in each candidate element.

Main inputs
	- Candidate regions BED file
	- DNase-seq of ATAC-seq input file
	- Genes reference file 

Output
	- EnhancerList.txt: Candidate enhancer regions with Dnase-seq (or ATAC-seq) and H3K27ac ChIP-seq read counts
	- GeneList.txt: Dnase-seq (or ATAC-seq) and H3K27ac ChIP-seq read counts on gene bodies and gene promoter regions

Description: 
	- Counts DNase-seq (or ATAC-seq) and H3K27ac ChIP-seq reads in candidate enhancer regions

**Example:**

.. code-block:: console

	$ python workflow/scripts/run.neighborhoods.py \
		--candidate_enhancer_regions results/K562_chr22/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed  \
		--DHS example_chr/chr22/ENCFF860XAE.chr22.sorted.se.bam \
		--default_accessibility_feature DHS \
		--chrom_sizes reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv \
		--outdir results/K562_chr22/Neighborhoods \
		--genes example_chr/chr22/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.hg38.bed.sorted.uniq \
		--ubiquitously_expressed_genes reference/UbiquitouslyExpressedGenes.txt \
		--qnorm reference/EnhancersQNormRef.K562.txt \
		--H3K27ac example_chr/chr22/ENCFF790GFL.chr22.sorted.se.bam

Key concepts to explain in this section:
- Activity is estimated by read count
- Quantile normalization of activity
- Which assays perform well for estimating activity



2.1. Activity scales with read counts [Jesse to add]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

2.2. Quantile normalization [Jesse to add]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

2.3. Using different combinations of assays to estimate enhancer activity [Andreas to add]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
e.g. note differences in perofrmanc eofr ATAC, DHS, H3K27ac, possibly add the ENCODE activity assay figure here


3. Estimating enhancer-promoter 3D contact (Jesse to edit)
------------------------------------------
Intro about concept of estimating enhancer-promoter 3D contact frequency by counting reads in Hi-C

We have different ways of estimating contact
	- Cell type specific Hi-C
	- Cell type average Hi-C
	- Power-law function of distance

Importance of Hi-C and how a lot is coming as a function of distance


Example code for each:

3.1. Cell-type average Hi-C data (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

3.2. Cell-type specific Hi-C data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

3.3. Power-law function of distance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

3.4. Pipeline to Download Hi-C data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You'll need to download HiC data to a local directory via juicer

.. code-block:: console

	$ python workflow/scripts/juicebox_dump.py  \
		--hic_file https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined_30.hic \
		--juicebox "java -jar juicer_tools.jar" \
		--outdir example_chr22/input_data/HiC/raw/ \
		--chromosomes 22

Powerlaw will be fit to the HiC dir if you use snakemake. If you wish to fit manually, you can run

.. code-block:: console

	$ python src/compute_powerlaw_fit_from_hic.py \
		--hic_dir example_chr22/input_data/HiC/raw/ \
		--hic_type juicebox \
		--hic_resolution 5000 \
		--outDir example_chr22/input_data/HiC/raw/powerlaw/ \



4. Making predictions with different combinations of input datasets  (Andreas to add performance comparison plots)
------------------------------------------------------------------------

Note: This code should really be moved elsewhere e.g. to a new section, like 'computing the ABC score'.  Need to explain somewhere the denominator of the ABC score

Main inputs
	- EnhancerList.txt
	- GeneList.txt
	- Powerlaw params (from fitting powerlaw to HiC data)
	- HiC data

Output
	- EnhancerPredictionsAllPutative.txt.gz: Scores for enhancer gene pairs

Description: 
	- Makes predictions following the Activity by Contact model
	- Utilizes HiC data for contact; otherwise, uses powerlaw

**Example:**

.. code-block:: console

	$ python workflow/scripts/predict.py  \
		--enhancers results/K562_chr22/Neighborhoods/EnhancerList.txt \
		--outdir results/K562_chr22/Predictions \
		--powerlaw_params_tsv results/HiC_Powerlaw/b08206e1/hic.powerlaw.tsv \
		--score_column ABC.Score \
		--chrom_sizes reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv \
		--accessibility_feature DHS \
		--cellType K562_chr22 \
		--genes results/K562_chr22/Neighborhoods/GeneList.txt \
		--hic_dir example_chr/HiC_K562 \
		--hic_type juicebox \
		--hic_resolution 5000 \
		--scale_hic_using_powerlaw			                                                                                                            



5. Interpreting the ABC score (Andreas to add)
------------------------------------

- Benchmark against the CRISPR data
- Correlates with effect size, but not in a linear way
- Appropriate threshold are different for models that use different combinations of input datasets, and provided [here]

