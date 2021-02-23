#! /bin/bash

ml R/4.0.2-foss-2019b

for SEED in 1..10
do 
	sbatch --wrap="code/fit_cv_superlearner.R $SEED $1"
done