#! bin/bash 

input_file=$1
outdir=$2
codedir=$3
hicdir=$4
hic_resolution=$5
scale_hic_using_powerlaw=$6
threshold=$7
make_all_putative=$8

while read p;
do
	a=($p)
	echo ${a[0]} ${a[1]}
	python $codedir/src/predict.py --enhancers $outdir/Neighborhood_${a[0]}/EnhancerList.txt --genes $outdir/Neighborhood_${a[0]}/GeneList.txt --HiCdir $hicdir --hic_resolution $hic_resolution --scale_hic_using_powerlaw $scale_hic_using_powerlaw --threshold $threshold --cellType ${a[0]} --outdir $outdir/Prediction_${a[0]} --make_all_putative $make_all_putative 

done < $input_file

