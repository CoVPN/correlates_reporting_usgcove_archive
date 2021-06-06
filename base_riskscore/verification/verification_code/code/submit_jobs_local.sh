#! /bin/bash

for SEED in `seq 1 10`
do 
	Rscript code/fit_cv_superlearner.R $SEED $1 > log$SEED.Rout &
done