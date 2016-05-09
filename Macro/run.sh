#!/bin/bash

source /lustre/nyx/cbm/users/klochkov/cbmroot/24_03_2016/build/config.sh

for run_id in `cat /lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Macro/RunsList.txt`
do
    echo $run_id
    root -l -b -q "TreeInterface.C  ($run_id)" &> log_$run_id.txt  &
done
