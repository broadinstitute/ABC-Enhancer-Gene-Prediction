#!/bin/bash

CODEDIR=/oak/stanford/groups/engreitz/Users/atan5133/ABC-Enhancer-Gene-Prediction

snakemake -s $CODEDIR/workflow/Snakefile_working --configfile $CODEDIR/config/config.yaml --directory $CODEDIR/workflow --unlock 

snakemake -s $CODEDIR/workflow/Snakefile_working --configfile $CODEDIR/config/config.yaml --use-conda -j1 --conda-frontend conda --keep-target-files --rerun-incomplete