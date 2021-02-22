## clean_data           : prepare a data set for SuperLearner and save in data_clean/ph1_placebo_final.rds
clean_data: code/data_prep.R
	Rscript code/data_prep.R

## seeds                : random ten seeds to set for Cv-Superlearner runs saved in output/seeds.rds
seeds: code/get_cv_sl_seeds.R
	Rscript code/get_cv_sl_seeds.R

## cvsl_out_slurm_demo : submit jobs to cluster to run cv super learner with demo library
cvsl_out_slurm_demo: code/fit_cv_superlearner.R seeds clean_data
	chmod +x code/fit_cv_superlearner.R && chmod +x code/submit_jobs.sh && code/submit_jobs.sh demo

## cvsl_out_slurm_prod : submit jobs to cluster to run cv super learner with production library
cvsl_out_slurm_prod: code/fit_cv_superlearner.R seeds clean_data
	chmod +x code/fit_cv_superlearner.R && chmod +x code/submit_jobs.sh && code/submit_jobs.sh prod