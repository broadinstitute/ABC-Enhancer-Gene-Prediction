.. _ABC-methods:

Detailed Methods
================

This section describes in more detail the key steps involved in calculating the ABC score, and discusses considerations around implementation choices that may be of interest to users and developers of ABC.

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
	    --accessibility example_chr/chr22/ENCFF860XAE.chr22.sorted.se.bam \
	    --outDir results/K562_chr22/Peaks \
	    --chrom_sizes reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv \
	    --chrom_sizes_bed results/tmp/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed \
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
	    --candidate_enhancer_regions results/K562_chr22/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
	    --DHS example_chr/chr22/ENCFF860XAE.chr22.sorted.se.bam \
	    --default_accessibility_feature DHS \
	    --chrom_sizes reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv \
	    --chrom_sizes_bed results/tmp/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed \
	    --outdir results/K562_chr22/Neighborhoods \
	    --genes results/K562_chr22/processed_genes_file.bed \
	    --ubiquitously_expressed_genes reference/UbiquitouslyExpressedGenes.txt \
	    --qnorm reference/EnhancersQNormRef.K562.txt \
	    --H3K27ac example_chr/chr22/ENCFF790GFL.chr22.sorted.se.bam

2.1. Activity scales with read counts 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Enhancer activity in the ABC model is estimated by counting reads in peaks (from DNase-seq, H3K27ac ChIP-seq, etc.) in peaks. The quantitative signal in these assays in informative regarding the strength of enhancers, and the ABC model assumes that this relationship is linear.

2.2. Quantile normalization for Activity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Datasets such as DNase-seq, ATAC-seq, and H3K27ac ChIP-seq often have varying signal-to-noise ratios (e.g., % reads in peaks, TSS enrichment). This changes the performance and thresholds needed for ABC model. To account for this, we apply quantile normalization on input datasets to match a reference dataset. As reference, we currently use datasets in K562, because we have CRISPR data to benchmark the model in that system.

2.3. Using different combinations of assays to estimate enhancer activity [Andreas to add]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
e.g. note differences in perofrmanc eofr ATAC, DHS, H3K27ac, possibly add the ENCODE activity assay figure here


3. Estimating enhancer-promoter 3D contact
------------------------------------------
The ABC model assumes that enhancer effects on gene expression vary linearly with 3D enhancer-promoter contact. In the initial version of the ABC model, we used quantitative observed signals from Hi-C datasets to estimate 3D contact (normalized by sequencing depth, but not by genomic distance). 

We now recommend selecting from one of three strategies to estimate enhancer-promoter 3D contact:

- *Cell-type specific Hi-C data*. If you have high-resolution Hi-C data available (e.g., 2+ billion reads for a genome-wide map), then using this data provides best performance for the model (see more details below)
- *A power law function of genomic distance*. If Hi-C data is not available, the simplest option that performs well is to estimate 3D contacts using a power law function of genomic distance. This power-law relationship explains >70% of the variance in Hi-C data (in situ Hi-C, 5-Kb resolution), and is sufficient for good performance, especially for shorter-range enhancer-gene pairs. This option should also be used when applying ABC to non-human organisms.
- *Cell-type average Hi-C data*. Another option for human samples is to use a cell-type averaged Hi-C map, in which the value for a given enhancer-promoter pair represents the average across available cell types. This method captures the relationship with genomic distance plus cell-type invariant 3D features such as certain long-range CTCF-mediated loops or domain boundaries.


Example biosample_config.tsv for each type:

3.1. Cell-type specific Hi-C data (best)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: auto

   * - biosample
     - DHS
     - ATAC
     - H3K27ac
     - default_accessibility_feature
     - HiC_file
     - HiC_type
     - HiC_resolution
     - alt_TSS
     - alt_genes
   * - K562
     - file/to/K562.bam
     - 
     - 
     - DHS
     - https://www.encodeproject.org/files/ENCFF621AIY/@@download/ENCFF621AIY.hic
     - hic
     - 5000
     - 
     - 

3.2. Power-law function of distance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: auto

   * - biosample
     - DHS
     - ATAC
     - H3K27ac
     - default_accessibility_feature
     - HiC_file
     - HiC_type
     - HiC_resolution
     - alt_TSS
     - alt_genes
   * - K562
     - file/to/K562.bam
     - 
     - 
     - DHS
     - 
     - 
     -
     - 
     - 

3.3. Cell-type average Hi-C data 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the avg HiC file here: https://www.encodeproject.org/files/ENCFF134PUN/@@download/ENCFF134PUN.bed.gz

To be filled out: Extract the file into multiple directories

.. list-table::
   :header-rows: 1
   :widths: auto

   * - biosample
     - DHS
     - ATAC
     - H3K27ac
     - default_accessibility_feature
     - HiC_file
     - HiC_type
     - HiC_resolution
     - alt_TSS
     - alt_genes
   * - K562
     - file/to/K562.bam
     - 
     - 
     - DHS
     - /path/to/avg_hic_directory
     - avg
     - 5000
     - 
     - 



3.4. Detailed considerations regarding estimation of 3D contact from Hi-C data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Starting from raw read counts in a Hi-C matrix, we perform several processing steps to normalize the data and handle missing or sparse data.

Normalization:

- We use SCALE or KR normalization to normalize by coverage in rows and columns
- For rows and columns from ENCODE Hi-C experiments corresponding to SCALE normalization factors < 0.25, we did not use SCALE normalization (these typically correspond to 5-kb bins with very few reads). Instead, we linearly interpolated the Hi-C signal in these bins by calculating an expected value based on power-law fit
- Each diagonal entry of the Hi-C matrix was replaced by the maximum of its four neighboring entries. This step is taken because the diagonal of the Hi-C contact map corresponds to the measured contact frequency between a 5-kb region of the genome and itself. The signal in bins on the diagonal can include restriction fragments that self-ligate to form a circle, or adjacent fragments that re-ligate, which are not representative of contact frequency. Empirically, we observed that the Hi-C signal in the diagonal bin was not well correlated with either of its neighboring bins and was influenced by the number of restriction sites contained in the bin.

We then compute Contact for an element-gene pair by rescaling the data as follows: 

- We set the Contact of the element-gene pair to the Hi-C signal at the bin of this row corresponding to the midpoint of E. For element-gene pairs that do not have a corresponding contact value (i.e., NaN), we set contact to zero. 
- For distances greater than 5 kb, we added a small adjustment (pseudocount) based on the power law expected count at a given distance threshold (as predicted by the power-law relationship between contact frequency and genomic distance). The distance threshold is usually determined by the resolution of the Hi-C data used. Distances less than 5 kb were given the pseudocount computed at 5 kb. We compute the power law by utilizing the steps in Section 6.3. 
- We found that different Hi-C datasets have slightly different power-law parameters. To weight all cell types equally in generating an average Hi-C profile, we scale the Hi-C profile in a given cell type by the cell-type specific gamma parameter from the power law relationship in that cell type. The scaling factor at distance d is given by d ^ (gamma_ref – gamma_celltype), where gamma_ref is the reference gamma parameter. 

Other considerations:

- We currently use Hi-C data at 5-Kb resolution.  Note that increasing the resolution is expected to affect the performance of the model both because of the potential sparsity in the data and because of the approach to normalizing Contact close to the TSS described above.


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

	$ python workflow/scripts/predict.py \
	    --enhancers results/K562_chr22/Neighborhoods/EnhancerList.txt \
	    --outdir results/K562_chr22/Predictions \
	    --score_column ABC.Score \
	    --chrom_sizes reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv \
	    --accessibility_feature DHS \
	    --cellType K562_chr22 \
	    --genes results/K562_chr22/Neighborhoods/GeneList.txt \
	    --hic_gamma 1.024238616787792 \
	    --hic_scale 5.9594510043736655 \
	    --hic_file https://www.encodeproject.org/files/ENCFF621AIY/@@download/ENCFF621AIY.hic \
	    --hic_type hic \
	    --hic_resolution 5000 \
	    --scale_hic_using_powerlaw			                                                                                                            



5. Interpreting the ABC score (Andreas to add)
------------------------------------

- Benchmark against the CRISPR data
- Correlates with effect size, but not in a linear way
- Appropriate threshold are different for models that use different combinations of input datasets, and provided [here]

