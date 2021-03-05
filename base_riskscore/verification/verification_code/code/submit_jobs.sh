#! /bin/bash

for SEED in `seq 1 10`
do 
	sbatch --output "log${SEED}.Rout" --wrap="code/fit_cv_superlearner.R $SEED $1"
done
