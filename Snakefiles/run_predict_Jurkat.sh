#!/bin/bash

# for Mike 
# THP1, hg19
python /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/scripts/predict_P2_ATAC.py \
--enhancers /oak/stanford/groups/engreitz/Projects/ABC/ABC_outputs/ABC_out/Jurkat-Engreitz/Neighborhoods/EnhancerList.txt \
--genes /oak/stanford/groups/engreitz/Projects/ABC/ABC_outputs/ABC_out/Jurkat-Engreitz/Neighborhoods/GeneList.txt \
--HiCdir /oak/stanford/groups/engreitz/Users/sheth/Data/ABC/hg19_average_hic \
--chrom_sizes /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/reference/hg19/chr_sizes \
--promoter_class /oak/stanford/groups/engreitz/Users/sheth/ABC2/data/for_classifier_p2class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.txt \
--hic_resolution 5000 \
--hic_type avg_hg19 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType THP1 \
--outdir /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/2023_1012_Jurkat \
--include_self_promoter \
--make_all_putative
