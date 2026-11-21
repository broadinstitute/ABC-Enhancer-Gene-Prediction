#!/bin/bash

# ABC_path=/oak/stanford/groups/engreitz/Users/kmualim/02212022_hg38Preds_QNORM-hic/K562_ID_2644
# code_dir=/oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction
# out_dir=$code_dir/2023_0115_Predictions_v1

# python $code_dir/Snakefiles/workflow/scripts/predict.py \
# --enhancers $ABC_path/Neighborhoods/EnhancerList.txt \
# --genes $ABC_path/Neighborhoods/GeneList.txt \
# --HiCdir /oak/stanford/groups/engreitz/Users/kmualim/for_maya/K562_intact_Hi-C \
# --chrom_sizes $code_dir/reference/hg38/hg38.chrom.sizes \
# --hic_resolution 5000 \
# --scale_hic_using_powerlaw \
# --threshold .02 \
# --cellType K562 \
# --outdir $out_dir \
# --include_self_promoter

python /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/scripts/predict_P2.py \
--enhancers /oak/stanford/groups/engreitz/Users/kmualim/02212022_hg38Preds_QNORM-hic/K562_ID_2644/Neighborhoods/EnhancerList.txt \
--genes /oak/stanford/groups/engreitz/Users/kmualim/02212022_hg38Preds_QNORM-hic/K562_ID_2644/Neighborhoods/GeneList.txt \
--promoter_class /oak/stanford/groups/engreitz/Users/sheth/ABC2/data/for_classifier_p2class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.txt \
--HiCdir /oak/stanford/groups/engreitz/Users/kmualim/for_maya/K562_intact_Hi-C \
--chrom_sizes /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/reference/hg38/hg38.chrom.sizes \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType K562 \
--outdir /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/2023_0619_ABC2_full_excludeP2_v1 \
--include_self_promoter \
--make_all_putative


# python /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/scripts/predict.py \
# --enhancers /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/Snakefiles/testEnhancerList.txt \
# --genes /oak/stanford/groups/engreitz/Users/kmualim/02212022_hg38Preds_QNORM-hic/K562_ID_2644/Neighborhoods/GeneList.txt \
# --HiCdir /oak/stanford/groups/engreitz/Users/kmualim/for_maya/K562_intact_Hi-C \
# --chrom_sizes /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/reference/hg38/hg38.chrom.sizes \
# --hic_resolution 5000 \
# --scale_hic_using_powerlaw \
# --threshold .02 \
# --cellType K562 \
# --outdir /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/2023_0201_Predictions_test \
# --include_self_promoter

# python /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/scripts/predict_P2.py \
# --enhancers $OAK/Users/kmualim/for_maya/K562_ENCDO000AAD_ENCFF325RTP_DNase-only/EnhancerList.txt \
# --genes $OAK/Users/kmualim/for_maya/K562_ENCDO000AAD_ENCFF325RTP_DNase-only/GeneList.txt \
# --HiCdir /oak/stanford/groups/engreitz/Users/kmualim/for_maya/K562_intact_Hi-C \
# --chrom_sizes /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/reference/hg38/hg38.chrom.sizes \
# --promoter_class /oak/stanford/groups/engreitz/Users/sheth/ABC2/data/for_classifier_p2class_RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.txt \
# --hic_resolution 5000 \
# --scale_hic_using_powerlaw \
# --threshold .02 \
# --cellType K562 \
# --outdir /oak/stanford/groups/engreitz/Users/sheth/ABC-Enhancer-Gene-Prediction-synergy/ABC-Enhancer-Gene-Prediction/2023_0617_ABC2_excludeP2_v1 \
# --include_self_promoter \
# --make_all_putative



