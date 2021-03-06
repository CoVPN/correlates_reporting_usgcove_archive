## clean_data           : prepare a data set for SuperLearner and save in data_clean/ph1_placebo_final.rds
clean_data: code/data_prep.R
	Rscript code/data_prep.R

## seeds                : random ten seeds to set for Cv-Superlearner runs saved in output/seeds.rds
seeds: code/get_cv_sl_seeds.R
	Rscript code/get_cv_sl_seeds.R

## cvsl_out_slurm_demo : submit jobs to cluster to run cv super learner with demo library
cvsl_out_slurm_demo: code/fit_cv_superlearner.R seeds clean_data code/submit_jobs.sh
	chmod +x code/fit_cv_superlearner.R && chmod +x code/submit_jobs.sh && code/submit_jobs.sh demo

## cvsl_out_local_demo : run jobs locally to fit cv super learner with demo library
cvsl_out_local_demo: code/fit_cv_superlearner.R seeds clean_data code/submit_jobs_local.sh
	chmod +x code/fit_cv_superlearner.R && chmod +x code/submit_jobs_local.sh && code/submit_jobs_local.sh demo

## cvsl_out_slurm_prod : submit jobs to cluster to run cv super learner with production library
cvsl_out_slurm_prod: code/fit_cv_superlearner.R seeds clean_data
	chmod +x code/fit_cv_superlearner.R && chmod +x code/submit_jobs.sh && code/submit_jobs.sh prod

## sl_out_local_demo : run jobs locally to fit cv super learner with demo library
sl_out_local_demo: code/fit_superlearner.R seeds clean_data code/submit_jobs_local.sh
	Rscript code/fit_superlearner.R 1234 demo

## sl_out_slurm_prod : submit jobs to cluster to run cv super learner with production library
sl_out_slurm_prod: code/fit_superlearner.R seeds clean_data
	sbatch --output "sl_log${SEED}.Rout" --wrap="code/fit_superlearner.R 1234 prod"

## cvauc_table         :  summarize results of cv super learner and save table
cvauc_table: code/get_cvauc.R
	@echo "~~~~~~~~~~\nMake sure you've saved cvSuperLearner results\n~~~~~~~~~~\n" 
	Rscript code/get_cvauc.R

## report              : report results
report: report.Rmd
	Rscript -e "rmarkdown::render('report.Rmd')"

## clean               : remove contents of output and data_clean
clean:
	rm -f output/* data_clean/*