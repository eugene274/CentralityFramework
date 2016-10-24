#!/bin/bash

# Task name
#SBATCH -J Glauber

# Working directory on shared storage
#SBATCH -D /lustre/nyx/cbm/users/klochkov/soft/CentralityFramework/Glauber

# Standard and error output in different files
#SBATCH -o %j_%N.out.log
#SBATCH -e %j_%N.err.log
#SBATCH -t 0-08:00 # Runtime in D-HH:MM format

# Execute application code

source /lustre/nyx/cbm/users/klochkov/cbmroot/22_08_2016/build/config.sh

currentDir=`pwd`
echo "current dir:" $currentDir

root -l -b -q "RunGlauberFitter.C  ($1, $2, $3, $4, $5, $6, $7)" || exit 11


echo "Done!"

