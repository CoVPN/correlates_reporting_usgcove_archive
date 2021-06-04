## all                  : get analysis results, but do not build report
all: cvauc_table risk_scores

## clean_data           : prepare a data set for SuperLearner and save in data_clean/ph1_placebo_final.rds
clean_data: code/data_prep.R
	Rscript code/data_prep.R

## seeds                : random ten seeds to set for Cv-Superlearner runs saved in output/seeds.rds
seeds: code/get_cv_sl_seeds.R
	Rscript code/get_cv_sl_seeds.R

## cvsl_out : run jobs locally to fit cv super learner with production library
cvsl_out: code/fit_cv_superlearner.R seeds clean_data code/submit_jobs_local.sh
	Rscript code/fit_cv_superlearner.R 1 prod

## sl_out_local_demo : run jobs locally to fit cv super learner with production library
sl_out: code/fit_superlearner.R seeds clean_data code/submit_jobs_local.sh
	Rscript code/fit_superlearner.R prod

## cvauc_table         :  summarize results of cv super learner and save table
cvauc_table: code/get_cvauc.R cvsl_out
	@echo "~~~~~~~~~~\nMake sure you've saved cvSuperLearner results\n~~~~~~~~~~\n" 
	Rscript code/get_cvauc.R

## risk_scores         : risk scores on the logit scale from superlearner fit
risk_scores: code/get_risk_scores.R sl_out
	Rscript code/get_risk_scores.R
	
## report              : report results
report: report.Rmd cvauc_table risk_scores
	Rscript -e "rmarkdown::render('report.Rmd')"

## clean               : remove contents of output and data_clean
clean:
	rm -f output/* data_clean/*

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<