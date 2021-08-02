## all           : clean saved results and get results to populate report
all: clean run_cvsl_riskscore createRDAfiles_fromSLobjects tables_figures \
	constructSL_predict_on_vaccine get_SLweights_Modelpredictors \
	append_risk_score

## clean         : delete all SL results
clean:
	rm -f output/*
	rm -f figs/*

## run_cvsl_riskscore   : run CV-SL on placebo arm
run_cvsl_riskscore:
	Rscript code/run_cvsl_riskscore.R
	#sbatch code/submit_cluster_job.sh
	
## createRDAfiles_fromSLobjects   : create RDA files from SL objects
createRDAfiles_fromSLobjects:
	Rscript code/createRDAfiles_fromSLobjects.R
	#sbatch code/createRDAfiles_fromSLobjects.sh
	
## tables_figures   : create tables and figures
tables_figures:
	Rscript code/tables_figures.R
	#sbatch code/tables_figures.sh
	
## constructSL_predict_on_vaccine   : construct SL and predict risk score on vaccinees
constructSL_predict_on_vaccine:
	Rscript code/constructSL_predict_on_vaccine.R
	#sbatch code/constructSL_predict_on_vaccine.sh
	
## get_SLweights_Modelpredictors   : get SL weights from models
get_SLweights_Modelpredictors:
	Rscript code/get_SLweights_Modelpredictors.R
	#sbatch code/get_SLweights_Modelpredictors.sh

## append_risk_score         : appends the computed risk scores to cleaned data
append_risk_score:
	Rscript code/append_risk_score_to_data.R

## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<
