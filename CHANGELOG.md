# Changelog


## 1.0.0 (2023-12-05)

### New Features

Snakemake Pipeline
* All the different components of ABC were moved into Snakemake, allowing easy execution of the entire ABC codebase with 1 snakemake command

Customization of ABC via config files
* Customization of ABC done in config files, instead of CLI arguments

Streaming of .hic files
* Support reading .hic files directly for the contact portion of the model

Support ATAC and scATAC as input 

Support passing multiple input files (e.g 2 DHS BAM files)

Documentation overhaul to ReadTheDocs

QC Plots for all ABC runs

Updated hg38 reference files

### Reliability Improvements

Conda Environments
* Conda environment yml file that can get built successfully on linux and macosx
* To prevent environment breakages, the conda environment gets built everyday via CircleCI

End to End Tests
* Tests that run ABC end to end and verifies there are no correctness regressions

### Bug Fixes/Improvements

Fixes to powerlaw computation, NaN ABC scores, sex chromosome normalization, and more

