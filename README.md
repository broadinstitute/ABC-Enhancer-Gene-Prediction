Check out [Read The Docs](https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/) for the most up to date information

# Activity by Contact Model of Enhancer-Gene Specificity

The Activity-by-Contact (ABC) model predicts which enhancers regulate which genes on a cell type specific basis. This repository contains the code needed to run the ABC model as well as small sample data files, example commands, and some general tips and suggestions. We provide a description of the model in the documentation, see Fulco, Nasser et al [1] for a full description.

v1.0 is the  most updated version of the code, with improved usability and performance. Usage is a bit different than previous version, so make sure to read the [getting started](https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/usage/getting_started.html#) page.

To use the version described in [1], we recommend release v0.2. There are some minor methodological differences between v0.2 and the model as described in [1]. These differences are related to Hi-C data processing and were implemented to improve the speed and scalability of the codebase. As such ABC scores computing using v0.2 will not exactly match those published in [1], although they will be very close. The codebase used to generate the results in [1] is available in the NG2019 branch of this repository. The NG2019 branch is no longer maintained.

If you use the ABC model in published research, please cite:

[1] Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). https://www.nature.com/articles/s41588-019-0538-0

## Contact
Please submit a github issue with any questions or if you experience any issues/bugs. 


