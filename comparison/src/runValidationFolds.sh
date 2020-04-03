#! bin/bash 
while read p;
do
	a=($p)
	python GrabValidationFolds.py --preds /mnt/lab_data2/kmualim/modified_ABC/results/mar262020/Predictions_${a[0]}_${a[3]}_IDR_candidateReg/EnhancerPredictionsAllPutative.txt.gz --outdir /mnt/lab_data2/kmualim/ABC_Comparison/${a[0]}_${a[3]}_IDR_candidateReg --expt /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/comparison/data/ExperimentalData.Gasperini.FulcoNasser.191021.txt
	python GrabValidationFolds.py --preds /mnt/lab_data2/kmualim/modified_ABC/results/mar262020/Predictions_${a[0]}_${a[3]}_optimal_IDR_narrowPeak/EnhancerPredictionsAllPutative.txt.gz --outdir /mnt/lab_data2/kmualim/ABC_Comparison/${a[0]}_${a[3]}_optimal_IDR_narrowPeak --expt /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/comparison/data/ExperimentalData.Gasperini.FulcoNasser.191021.txt
done < $1 

