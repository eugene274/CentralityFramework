#!/bin/bash

firstfile=0
lastfile=2999
filesperjob=100

for filenum in `seq $firstfile $filesperjob $lastfile`
do
    echo $filenum 
    sbatch runCentralityConrainer.sh $filenum  $filenum+$filesperjob-1
done
