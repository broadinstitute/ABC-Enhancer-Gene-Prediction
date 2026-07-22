#!/bin/bash
#SBATCH --job-name=h3k27ac_counting_abc
#SBATCH -p normal,engreitz
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -o config/mine/logs/h3k27ac_counting_%j.out
#SBATCH -e config/mine/logs/h3k27ac_counting_%j.err

set -euo pipefail

cd /oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction

# Reuse the conda env snakemake already built for this repo (has snakemake,
# macs2, bedtools, samtools, pysam, pandas, pyranges, etc; see workflow/envs/abcenv.yml)
ENV=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/.snakemake/conda/a0febb704d007e3c12f577f28bdeccb8_
export PATH="$ENV/bin:$PATH"

# This driver job just orchestrates (low mem/cpu); --profile submits each
# rule as its own appropriately-sized sbatch job (see config/mine/slurm_profile).
snakemake -s workflow/Snakefile \
    --configfile config/mine/config_h3k27ac_counting.yaml \
    --profile config/mine/slurm_profile \
    -p
