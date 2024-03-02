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

The method of defining candidate elements includes the following steps:

- Peak-calling with MACS2
- Resizing and merging regions
- Selecting the 150,000 strongest peaks (by read count)
- Adding promoters

1.1. Calling peaks with MACS2 
------------------------------
Main inputs
	- ATAC-Seq:
		- Sorted tagalign file. The sorted tagalign file should be free of PCR duplicates, and sorted according to the chromosome order in the *chrom.sizes.tsv file; the file name must end in “.tagAlign.gz”. 
	- DNAse-Seq:
		- Single-ended BAM file. 

Output
	- narrowPeak file

Description:
	- To identify enriched cutting sites in DNase-Seq and ATAC-Seq datasets, we shift the 5’ ends of all reads by 75 bp towards the 3’ direction and extend them by 150 bp in the 5’ direction following the `ENCODE ATACSeq pipeline <https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit#heading=h.9ecc41kilcvq>`_.

	- Because (single-nucleus) ATAC-Seq is mostly paired-end, and MACS2 cannot call peaks using all the reads in paired-end BAM files, we recommend that ATAC-Seq input should be in the tagAlign format. Specifically, when using the “—format BAM” flag with paired-end BAM files, MACS only considers reads with R1 tags. Alternatively, the “—format BAMPE” flag would allow MACS2 to use all paired-end reads; however, “—format BAMPE”  cannot be performed with the reads shift and extension described above. 

	- Since most DNase-Seq is single ended, and the ”—format BAM” option is compatible with shift and extension, DNase-Seq input can be supplied in the BAM format.

1.2. Resizing and merging regions 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We resize to 500 bp to count reads in and around peaks (esp. e.g. H3K27ac signal is surrounding the peak)


1.3. Selecting the top N peaks
------------------------------
Description: 
	- To define the candidate regions, for genome-wide analyses, we retain the top 150,000 peaks with the most read counts. A fixed number is chosen here because the numbers of peaks called vary with sequencing depths, but imprically we discovered that picking the peaks with the most reads counts can effectively remove the noise coming from weak peaks and variable sequencing quality. Additionally, the number of total peaks also affect the denominator of ABC score calculation; a fixed number of peaks also make ABC scores comparable across inputs of variable sequencing qualities and depths. For genome-wide analyses, 150K is a reasonable number because ENCODE analysis has previously estimated `a mean of 205,109 DHSs per cell type <https://www.nature.com/articles/nature11247>`, the majority of which are enhancers. 

1.4. Defining and adding gene promoters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, we force the inclusion of gene promoters in the set of candidate elements, to include
promoters of all candidate genes the calculation of ABC scores. Promoters are chromatin accessible
sites, however sometimes the promoters of genes do not pass the threshold for top 150,000 strongest
peaks. Inclusion of promoter regions has large effect on the model due to the promoter receiving a
high "3D contact" value in the ABC computation.

Note that the exact method of defining the promoter region for a gene has a strong influence on
computing ABC and changing the promoter list can impact ABC scores:

#. The exact promoter is used to compute 3D contact for enhancers regulating a gene, and using a wrong promoter will result in incorrect 3D contact estimates. If the promoter that is included itself as a candidate element is the dominant regulatory elements for a gene, it heavily contributes to the denominator of the ABC score (see below). Therefore, using an incorrect promoter will impact the ABC scores of enhancers regulating a gene.
#. The list of promoters affects the ABC scores for other nearby genes, because force inclusion of these regions leads to more/larger regions being used as candidate elements, which affects the ABC denominator.
#. The exact promoter list can affect benchmarking and downstream analyses. For example, analyses that filter to include only 'distal elements' that are not promoters might filter out elements incorrectly called as promoters that are in reality enhancers (e.g. promoters of lncRNAs that act as enhancers).
	
In practice, we provide a gene promoter file that we have used for various purposes that selects a
single canonical promoter per gene. This file contains canonical RefSeq promoters which could be
assigned to a gene in the GENCODE v29 genome annotations used by the ENCODE consortium. Genes that
are annotated in GENCODE as miRNAs, pseudogenes, antisense transcribed RNAs or RNAs transcribed from
bidirectional promoters were filtered out.

Changing the promoter for a single gene, e.g. to accommodate a specific alternative transcription
start site of a gene of interest, is likely to not affect predictions globally and can be used in
certain cases. However, caution is warranted if making more extensive changes to the promoter list.

Also including a much larger promoter list, e.g. including lncRNAs or including all possible
transcription start sites for all isoforms for a gene, is likely to change the global properties of
the ABC score and is not recommended without calibration of scores (see section on Interpreting the
ABC score below).

Note that ABC includes the promoter of each gene in the thresholded ABC enhancer-gene predictions
regardless of it's ABC score (forced to 1), since the promoter of a gene is always considered to
regulate it's expression.



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

2.1. Activity scales with read counts 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Enhancer activity in the ABC model is estimated by counting reads in peaks (from DNase-seq, H3K27ac ChIP-seq, etc.) in peaks. The quantitative signal in these assays in informative regarding the strength of enhancers, and the ABC model assumes that this relationship is linear.

2.2. Quantile normalization for Activity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Datasets such as DNase-seq, ATAC-seq, and H3K27ac ChIP-seq often have varying signal-to-noise ratios (e.g., % reads in peaks, TSS enrichment). This changes the performance and thresholds needed for ABC model. To account for this, we apply quantile normalization on input datasets to match a reference dataset. As reference, we currently use datasets in K562, because we have CRISPR data to benchmark the model in that system.

2.3. Using different chromatin assays to estimate enhancer activity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ABC uses DNase-seq and optionally H3K27ac ChIP-seq to estimate enhancer activity, but numerous other
chromatin assays exist, including chromatin accessibility assays, TF ChIP-seq and histone ChIP-seq.
To investigate which assays perform best in the ABC framework, we've built ABC models where activity
was measured by one of 513 ENCODE 1D chromatin experiments that could represent enhancer activity.

Among all assays, ABC models using DNase-seq or H3K27ac ChIP-seq ranked among the top 10 when
benchmarking the performance of these models against CRISPR enhancer perturbation results. Given the
availability of DNase-seq and H3K27ac ChIP-seq datasets (e.g. www.encodeproject.org), these two
assays provide the best performance to build models that can be applied across different cell types
and tissues. Of note, bulk ATAC-seq performed worse than DNase-seq in this comparison.

.. image:: /images/abc_enh_assays_perf.png


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

We need to extract the file to create multiple HiC directories, one for each chromosome. 

.. code-block:: console

  (abc-env) [atan5133@sh03-04n24 /oak/stanford/groups/engreitz/Users/atan5133/ABC-Enhancer-Gene-Prediction] (job 37343981) $ py workflow/scripts/extract_avg_hic.py --avg_hic_bed_file ../data/ENCFF134PUN.bed.gz --output_dir ../data/

This may take a while (~1-2 hours)

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
     - /oak/stanford/groups/engreitz/Users/atan5133/data/AvgHiC
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


4. Making predictions with different combinations of input datasets
------------------------------------------------------------------------

ABC scores are computed for all generated candidate elements within 5Mb of the considered
candidate target genes (GeneList.txt). For each element-pair, the activity of the element is
multiplied by the 3D contact value between the element and the promoter. 3D contact values can be
calculated using different approaches and we achieved the highest performance by using cell-type
specific Hi-C data. If Hi-C isn't available for the cell type of interest, average Hi-C contact
across cell types (described in
`Nasser et al., 2021 <https://www.nature.com/articles/s41586-021-03446-x>`_) or an approximation of
3D contact via a powerlaw function (described in
`Fulco et al., 2019 <https://www.nature.com/articles/s41588-019-0538-0>`_) can be used to calculate
contact.

Finally, the ABC scores for all element-gene pairs are divided by the sum of all ABC scores for each
gene, normalizing the sum of ABC scores per gene to 1.

The precision-recall curves below show a comparison of the CRISPR benchmark performance of ABC
models using cell-type specific Hi-C versus the powerlaw approximation for elements inferred from
bulk DNase-seq and single-cell ATAC-seq:

.. image:: /images/abc_perf_comparison.png


Main inputs
	- EnhancerList.txt
	- GeneList.txt
	- Powerlaw params (from fitting powerlaw to HiC data)
	- HiC data (optional)

Output
	- EnhancerPredictionsAllPutative.txt.gz: Scores for enhancer gene pairs

Description: 
	- Makes predictions following the Activity by Contact model
	- Utilizes HiC data for contact; otherwise, uses powerlaw

5. Interpreting the ABC score
------------------------------------

To validate and better understand the properties of the ABC score, we extensively benchmarked the
model against CRISPR enhancer perturbation in K562 cells (
`Fulco et al., 2019 <https://www.nature.com/articles/s41588-019-0538-0>`_ ,
`Nasser et al., 2021 <https://www.nature.com/articles/s41586-021-03446-x>`_).

These analyses show that ABC scores reliably predicts enhancer-gene regulatory interactions
that were experimentally inferred in the CRISPR experiments. At the recall of 70%, an ABC model
using DNase-seq + cell-type specific Hi-C data achieves a precision of 51%, meaning around half of
the predicted enhancer-gene regulatory interactions will be true positives. The ABC scores
themselves correlate with the CRISPR effect size on gene expression when perturbing an enhancer,
however not in a precise linear fashion. This probably has different technical and biological
reasons. Nevertheless, we can expect enhancers with large ABC scores to have strong effects on gene
expression.

One key consideration when applying ABC to generate maps of enhancer-gene pairs is selecting the
appropriate ABC score threshold to predict regulatory enhancer-gene interactions. We recommend using
a threshold that achieves 70% recall in our CRISPR benchmark. A list of thresholds for different
common ABC models can be found in the table below:

.. csv-table::
   :file: /tables/perf_summary.csv
   :header-rows: 1
    
We automatically choose the best threshold based on your input, but you can specify a threshold value
yourself in the config.yaml file.

Our CRISPR benchmarking pipeline can be used to infer thresholds for non-standard ABC models and is
available on `Github <https://github.com/EngreitzLab/CRISPR_comparison>`_.
