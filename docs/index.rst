.. ABC-Enhancer-Gene-Prediction documentation master file, created by
   sphinx-quickstart on Mon Jun 12 16:07:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Activity-by-Contact (ABC) Model for predicting enhancer-gene regulatory interactions
========================================================================================

The Activity-by-Contact (ABC) model predicts which enhancers regulate which genes in a given cell type, based on a combination of input datasets representing enhancer activity and 3D enhancer-promoter contact frequency [1].

The ABC model was initially designed to represent a mechanistic model in which enhancers activate gene transcription upon enhancer-promoter contact. In a simple conception of such a model, the quantitative effect of an enhancer depends on the frequency with which it contacts a promoter multiplied by the strength of the enhancer (i.e., the ability of the enhancer to activate transcription upon contacting a promoter). Moreover, the contribution of a specific enhancer to a gene’s expression should depend on the surrounding context (ie, the strength and contact frequency of other enhancers for the gene).  To convert this conceptual framework into a practical score (which can be applied genome-wide), we formulated the ABC score:

ABC score for effect of element E on gene G = Activity of E × Contact frequency between E and G / Sum of (Activity × Contact Frequency) over all candidate elements within 5 Mb.

Operationally, Activity (A) is defined based on read counts of ATAC-seq, DNase-seq, and/or H3K27ac ChIP-seq at an element E, and Contact (C) as the Hi-C contact frequency between E and the promoter of gene G. Elements are defined as ~500bp regions centered on DHS or ATAC peaks.

Note that the ABC model only considers candidate elements and genes on the same chromosome, within 5 Mb of each other. It does not make interchromosomal predictions.

The accuracy of the ABC model has been benchmarked in several ways, including through comparison to CRISPR perturbations to enhancers, eQTL variants, and GWAS variants. For more information, please see the references below.

If you use the ABC model in published research, please cite:

[1] Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). https://www.nature.com/articles/s41588-019-0538-0

[2] Nasser J, Bergman DT, Fulco CP, Guckelberger P, Doughty BR, Patwardhan TA, Jones TR, Nguyen TH, Ulirsch JC, Lekschas F, Mualim K, Natri HM, Weeks EM, Munson G, Kane M, Kang HY, Cui A, Ray JP, Eisenhaure TM, Collins RL, Dey K, Pfister H, Price AL, Epstein CB, Kundaje A, Xavier RJ, Daly MJ, Huang H, Finucane HK, Hacohen N, Lander ES, Engreitz JM. Genome-wide enhancer maps link risk variants to disease genes. Nature. 2021 May;593(7858):238-243. doi: 10.1038/s41586-021-03446-x

Released Versions
-----------------

v1.0 (version for ENCODE paper)

- Snakemake and config files for improved usability
- QC Plots
- Other feature additions and bug fixes

v0.2 (version for [1])

- There are some minor methodological differences between v0.2 and the model as described in [1]. These differences are related to Hi-C data processing and were implemented to improve the speed and scalability of the codebase. As such ABC scores computing using v0.2 will not exactly match those published in [1], although they will be very close. The codebase used to generate the results in [1] is available in the NG2019 branch of this repository. The NG2019 branch is no longer maintained.



.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Usage

   usage/getting_started
   usage/methods
   usage/scATAC
   usage/troubleshooting


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Developers

   developers/contributing
   developers/testing
   developers/updating_docs
   developers/internal_docs

