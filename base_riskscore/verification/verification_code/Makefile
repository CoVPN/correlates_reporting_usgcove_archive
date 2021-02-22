## clean_data     : prepare a data set for SuperLearner and save in data_clean/ph1_placebo_final.rds
clean_data: code/data_prep.R
	Rscript code/data_prep.R

## seeds          : random ten seeds to set for Cv-Superlearner runs saved in output/seeds.rds
seeds: code/get_cv_sl_seeds.R
	Rscript code/get_cv_sl_seeds.R

