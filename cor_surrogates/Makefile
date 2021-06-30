## all           : clean saved results and get results to populate report
all: clean run_cvsl_varsets createRDAfiles_fromSLobjects tables_figures

## clean         : delete all SL results
clean:
	rm -f output/*
	rm -f results/*

## run_cvsl_varsets   : run CV-SL on vaccine arm and each variable set
run_cvsl_varsets: 
	Rscript code/run_cvsl_varsets.R 1 && \
	Rscript code/run_cvsl_varsets.R 2 && \
	Rscript code/run_cvsl_varsets.R 3 && \
	Rscript code/run_cvsl_varsets.R 4 && \
	Rscript code/run_cvsl_varsets.R 5 && \
	Rscript code/run_cvsl_varsets.R 6 && \
	Rscript code/run_cvsl_varsets.R 7 && \
	Rscript code/run_cvsl_varsets.R 8 && \
	Rscript code/run_cvsl_varsets.R 9 && \
	Rscript code/run_cvsl_varsets.R 10 && \
	Rscript code/run_cvsl_varsets.R 11 && \
	Rscript code/run_cvsl_varsets.R 12
	#sbatch code/submit_cluster_job.sh
	
## createRDAfiles_fromSLobjects   : create RDA files from SL objects
createRDAfiles_fromSLobjects: run_cvsl_varsets
	Rscript code/createRDAfiles_fromSLobjects.R
	#sbatch code/createRDAfiles_fromSLobjects.sh
	
## tables_figures   : create tables and figures
tables_figures: createRDAfiles_fromSLobjects
	Rscript code/tables_figures.R
	#sbatch code/tables_figures.sh
	
## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<
