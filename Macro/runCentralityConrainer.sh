#!/bin/bash

# Task name
#SBATCH -J CentrCont

# Working directory on shared storage
#SBATCH -D /lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Macro

# Standard and error output in different files
#SBATCH -o %j_%N.out.log
#SBATCH -e %j_%N.err.log
#SBATCH -t 0-01:00 # Runtime in D-HH:MM format

# Execute application code

source /lustre/nyx/cbm/users/klochkov/cbmroot/24_03_2016/build/config.sh

currentDir=`pwd`
echo "current dir:" $currentDir

root -l -b -q "TreeInterface.C  ($1)"

echo "Done!"

