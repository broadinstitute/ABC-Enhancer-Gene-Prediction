#!/bin/sh
# properties = {properties}

source /broad/software/scripts/useuse
reuse -q .python-3.5.1
source /seq/lincRNA/Ben/VENV_BEN_UMITOOLS/bin/activate
reuse -q BEDTools
reuse -q Bowtie
reuse -q Samtools
reuse -q FASTX-Toolkit
reuse BLAST
# reuse Matlab
reuse .matlab-2017a

{exec_job}
