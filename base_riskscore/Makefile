## all           : clean saved results and get results to populate report
all: clean run_cvsl_riskscore createRDAfiles_fromSLobjects tables_figures constructSL_predict_on_vaccine get_SLweights_Modelpredictors


## clean         : delete all SL results
clean:
	rm -f output/*
	rm -f results/*

## run_cvsl_riskscore   : run CV-SL on placebo arm
run_cvsl_riskscore:
	sbatch code/submit_cluster_job.sh 
	
## createRDAfiles_fromSLobjects   : create RDA files from SL objects
createRDAfiles_fromSLobjects:
	sbatch code/createRDAfiles_fromSLobjects.sh 
	
## tables_figures   : create tables and figures
tables_figures:
	sbatch code/tables_figures.sh 
	
## constructSL_predict_on_vaccine   : construct SL and predict risk score on vaccinees
constructSL_predict_on_vaccine:
	sbatch code/constructSL_predict_on_vaccine.sh 
	
## get_SLweights_Modelpredictors   : get SL weights from models
get_SLweights_Modelpredictors:
	sbatch code/get_SLweights_Modelpredictors.sh 

## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_coxph
