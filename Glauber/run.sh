#!/bin/bash

nf=1
nsigma=15
step_sigma=0.1
n_events=2000000
n_f_steps=50
min_mult=$1

mkdir root_files/MinMult_$min_mult

for i in `seq 0 $n_f_steps`;
do
#     echo $i
    f=1.0*$i/$n_f_steps
    echo $f
    sbatch run_glauber.sh $nf $f $f $nsigma $step_sigma $n_events $min_mult
done 
