#!/bin/bash
# Job name:
#SBATCH --job-name=covpn_cop_stochastic
#
# Working directory:
#SBATCH --chdir=/global/home/users/nhejazi/
#
# Account:
#SBATCH --account=co_biostat
#
# Partition:
#SBATCH --partition=savio3
#
# Processors (1 node = 20 cores):
#SBATCH --nodes=1
#SBATCH --exclusive
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=20:00:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=nhejazi@berkeley.edu
#
# Job output:
#SBATCH --output=slurm.out
#SBATCH --error=slurm.out
#
## Command(s) to run:
export TMPDIR='/global/scratch/nhejazi/rtmp'  # resolve update issues for compiled packages as per https://github.com/r-lib/devtools/issues/32
export R_LIBS_USER='/global/scratch/nhejazi/R'  # personal package library
module load gcc/6.3.0 r/4.0.3 r-packages/default
cd ~/correlates_reporting
make data_processed
echo "run_fast <- FALSE" >> ./cop_stochastic/code/params.R
make -k -C cop_stochastic all
