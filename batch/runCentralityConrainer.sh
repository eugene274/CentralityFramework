#!/bin/bash

# Task name
#SBATCH -J CentrCont

# Working directory on shared storage
#SBATCH -D /lustre/nyx/cbm/users/klochkov/git/CentralityFramework/batch

# Standard and error output in different files
#SBATCH -o %j_%N.out.log
#SBATCH -e %j_%N.err.log
#SBATCH -t 0-01:00 # Runtime in D-HH:MM format

# Execute application code

source /lustre/nyx/cbm/users/klochkov/cbm/mc/DCM_QGSM/Au4Au/sts_psd_only/20161221/cbmroot/build/config.sh

currentDir=`pwd`
echo "current dir:" $currentDir

root -l -b -q "/lustre/nyx/cbm/users/klochkov/git/CentralityFramework/Macro/RunTreeInterface.C  ($1, $2)"

echo "Done!"

