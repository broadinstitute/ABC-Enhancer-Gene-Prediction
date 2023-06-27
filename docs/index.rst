.. ABC-Enhancer-Gene-Prediction documentation master file, created by
   sphinx-quickstart on Mon Jun 12 16:07:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ABC-Enhancer-Gene-Prediction's documentation!
========================================================

The Activity-by-Contact (ABC) model predicts which enhancers regulate which genes on a cell type specific basis.

If you use the ABC model in published research, please cite:

[1] Fulco CP, Nasser J, Jones TR, Munson G, Bergman DT, Subramanian V, Grossman SR, Anyoha R, Doughty BR, Patwardhan TA, Nguyen TH, Kane M, Perez EM, Durand NC, Lareau CA, Stamenova EK, Aiden EL, Lander ES & Engreitz JM. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019). https://www.nature.com/articles/s41588-019-0538-0

[2] Nasser J, Bergman DT, Fulco CP, Guckelberger P, Doughty BR, Patwardhan TA, Jones TR, Nguyen TH, Ulirsch JC, Lekschas F, Mualim K, Natri HM, Weeks EM, Munson G, Kane M, Kang HY, Cui A, Ray JP, Eisenhaure TM, Collins RL, Dey K, Pfister H, Price AL, Epstein CB, Kundaje A, Xavier RJ, Daly MJ, Huang H, Finucane HK, Hacohen N, Lander ES, Engreitz JM. Genome-wide enhancer maps link risk variants to disease genes. Nature. 2021 May;593(7858):238-243. doi: 10.1038/s41586-021-03446-x

Released Versions (To be updated)
----------------------------------

v0.2 is the recommended version for the majority of users. There are some minor methodological differences between v0.2 and the model as described in [1]. These differences are related to Hi-C data processing and were implemented to improve the speed and scalability of the codebase. As such ABC scores computing using v0.2 will not exactly match those published in [1], although they will be very close. The codebase used to generate the results in [1] is available in the NG2019 branch of this repository. The NG2019 branch is no longer maintained.



.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Usage

   usage/getting_started
   usage/about_abc


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Developers

   developers/developers

