Detailed Methods
=========

To-do:  This page should include detailed methods for each of the main operations for ABC

We may want to include benchmarking figures or other methodological details e.g. from previous papers

Key concepts:

- Defining candidate elements
- Estimating enhancer activity
- Estimating enhancer-promoter 3D contact
- Making predictions with different combinations of input datasets
- Interpreting the ABC score

1. Defining candidate elements
---------

'Candidate elements' are the set of putative enhancer elements for which ABC Scores will be computed. These include gene promoters. 

The method of defining candidate elements includes the following steps:

- Peak-calling with MACS2
- Resizing and merging regions
- Selecting the 150,000 strongest peaks (by read count)
- Adding promoters

1.1. Calling peaks with MACS2
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Rosa to add details about the MACS2 methods calling

  
1.2. Defining gene promoters
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Defining the promoter region for a gene has a strong influence on the ABC computation.

describe how it is important how promoters are selected, and how changing the promoter list can impact ABC scores 
- First, the exact promoter used affects the ABC score for the gene corresponding to that promoter, because of 3D contacts (which can differ depending on the location of the promoter) and whether that promoter is in fact the dominant element used (the promoter is included as a candidate "enhancer" for itself, and contributes to the denominator of the ABC score)
- Second, the promoter list used affects the ABC scores for other nearby genes, because force inclusion of these regions leads to more/larger regions being used which affects the ABC denominator
- Third, the promoter list used can affect downstream benchmarking analyses. For example, benchmarks that filter to just 'distal elements' that are not promoter might filter out elements called as promoters that are actually enhancers (e.g. promoters of lncRNAs that act as enhancers)

In practice, we provide a gene promoter file that we have used for various purposes that selects a single canonical promoter per gene. 
- describe provenance of the gene promoter file(s) including in ABC repo (for human and mouse)
- Changing the promoter for a single gene, e.g. to accommodate a specific alternative transcription start site of a gene of interest, is likely to be okay and not globally affect predictions
- However, caution is warranting in making more extensive changes to the promoter list. Note again that including a much larger promoter list, e.g. including lncRNAs or including all possible transcription start sites for all isoforms for a gene, is likely to change the global properties of the ABC score and is not recommended without calibration of scores (see section on Interpreting the ABC score below)


2. Estimating enhancer activity
---------

3. Estimating enhancer-promoter 3D contact
---------

3.1. Cell-type average Hi-C data (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

3.2. Cell-type specific Hi-C data
^^^^^^^^^^^^^^^^^^^^^^^^^^^

3.3. Power-law function of distance
^^^^^^^^^^^^^^^^^^^^^^^^^^^


4. Making predictions with different combinations of input datasets
---------


5. Interpreting the ABC score
---------

- Benchmark against the CRISPR data
- Correlates with effect size, but not in a linear way
- Appropriate threshold are different for models that use different combinations of input datasets, and provided [here]

