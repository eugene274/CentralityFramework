#!/bin/bash
run_id=23005
# for run_id in `cat /lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Macro/RunsList.txt`
# do
    echo $run_id
    sbatch runCentralityConrainer.sh $run_id
# done
