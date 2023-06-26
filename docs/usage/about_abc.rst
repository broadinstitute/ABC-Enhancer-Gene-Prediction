About ABC
=========

The Activity by Contact (ABC) model is designed to represent a mechanistic model in which enhancers activate gene transcription upon enhancer-promoter contact. In a simple conception of such a model, the quantitative effect of an enhancer depends on the frequency with which it contacts a promoter multiplied by the strength of the enhancer (i.e., the ability of the enhancer to activate transcription upon contacting a promoter). Moreover, the contribution of a specific enhancer to a gene’s expression should depend on the surrounding context (ie, the strength and contact frequency of other enhancers for the gene).

To convert this conceptual framework into a practical score (which can be applied genome-wide), we formulated the ABC score:

ABC score for effect of element E on gene G = Activity of E × Contact frequency between E and G / Sum of (Activity × Contact Frequency) over all candidate elements within 5 Mb.

Operationally, Activity (A) is defined as the geometric mean of the read counts of DNase-seq and H3K27ac ChIP-seq at an element E, and Contact (C) as the KR normalized Hi-C contact frequency between E and the promoter of gene G. Elements are defined as ~500bp regions centered on DHS peaks.

Note that the ABC model only considers candidate elements and genes on the same chromosome. It does not make interchromosomal predictions.

.. _ABC-steps:

ABC Steps
---------

These are the high level steps performed when running the ABC pipeline

Step 1. Define candidate elemets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

'Candidate elements' are the set of putative enhancer elements for which ABC Scores will be computed. A typical way to define candidate elements is by calling peaks on a DNase-Seq or ATAC-Seq bam file. In this implementation we first call peaks using MACS2 and then process these peaks using ``makeCandidateRegions.py``.

``makeCandidateRegions.py`` will take as input the narrowPeak file produced by MACS2 and then perform the following processing steps:

#. Count DNase-seq reads in each peak and retain the top N peaks with the most read counts
#. Resize each of these N peaks to be a fixed number of base pairs centered on the peak summit
#. Remove any regions listed in the 'blocklist' and include any regions listed in the 'includelist'
#. Merge any overlapping regions

Main output files:

- *macs2_peaks.narrowPeak*: MACS2 narrowPeak file
- *macs2_peaks.narrowPeak.candidateRegions.bed*: filtered, extended and merged peak calls from MACS2. These are the candidate regions used in downstream scripts.

Step 2. Quantifying Enhancer Activity:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``run.neighborhoods.py`` will count DNase-seq (or ATAC-seq) and H3K27ac ChIP-seq reads in candidate enhancer regions. It also makes GeneList.txt, which counts reads in gene bodies and promoter regions.

Replicate epigenetic experiments should be included as comma delimited list of files. Read counts in replicate experiments will be averaged when computing enhancer Activity.

Main output files:

  - **EnhancerList.txt**: Candidate enhancer regions with Dnase-seq and H3K27ac ChIP-seq read counts
  - **GeneList.txt**: Dnase-seq and H3K27ac ChIP-seq read counts on gene bodies and gene promoter regions

Step 3. Computing the ABC Score
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute ABC scores by combining Activity (as calculated by ``run.neighborhoods.py``) and Hi-C.


