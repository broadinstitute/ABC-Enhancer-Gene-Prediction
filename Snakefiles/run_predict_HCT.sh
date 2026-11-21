#!/bin/bash

# for Philline, 
# NEighborhoods are the same (only Hi-C differs btw treated and untreated)
# UNTREATED
python /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/scripts/predict.py \
--enhancers /oak/stanford/groups/engreitz/Users/kmualim/for_maya/HCT116_ENCFF166PCV_ENCFF799ZUN_ENCFF977FPK_ENCFF938MXK/Neighborhoods/EnhancerList.txt \
--genes /oak/stanford/groups/engreitz/Users/kmualim/for_maya/HCT116_ENCFF166PCV_ENCFF799ZUN_ENCFF977FPK_ENCFF938MXK/Neighborhoods/GeneList.txt \
--HiCdir /oak/stanford/groups/engreitz/Projects/ABC/HiC/HCT116-RAD21_auxin_untreated_intact_HIC \
--chrom_sizes /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/reference/hg38/hg38.chrom.sizes \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType HCT116-untreated \
--outdir /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/2023_0518_HCT116_untreated \
--include_self_promoter \
--make_all_putative

# TREATED
python /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/scripts/predict.py \
--enhancers /oak/stanford/groups/engreitz/Users/kmualim/for_maya/HCT116_ENCFF166PCV_ENCFF799ZUN_ENCFF977FPK_ENCFF938MXK/Neighborhoods/EnhancerList.txt \
--genes /oak/stanford/groups/engreitz/Users/kmualim/for_maya/HCT116_ENCFF166PCV_ENCFF799ZUN_ENCFF977FPK_ENCFF938MXK/Neighborhoods/GeneList.txt \
--HiCdir /oak/stanford/groups/engreitz/Projects/ABC/HiC/HCT116-RAD21-treated \
--chrom_sizes /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/reference/hg38/hg38.chrom.sizes \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType HCT116-untreated \
--outdir /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/2023_0518_HCT116_untreated \
--include_self_promoter \
--make_all_putative