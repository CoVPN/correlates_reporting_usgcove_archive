#! /bin/bash

for SEED in 1..10
do 
	sbatch --wrap="code/fit_cv_superlearner.R $SEED $1"
done