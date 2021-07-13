## all           : clean saved results and get results to populate report
all: clean run_cvsl_varsets_d29 createRDAfiles_fromSLobjects_d29 tables_figures_d29 \
	run_cvsl_varsets_d57 createRDAfiles_fromSLobjects_d57 tables_figures_d57 \
	run_cvsl_varsets_Both createRDAfiles_fromSLobjects_Both tables_figures_Both 

## clean         : delete all SL results
clean:
	rm -f output/*
	rm -f figs/*

######################################################################
## DAY 29 MARKERS ANALYSES
run_cvsl_varsets_d29:
	for number in 1 2 3 4 5 6 7 8 9 10; do \
	Rscript code/run_cvsl_varsets.R $$number Day29; \
	done
    
## createRDAfiles_fromSLobjects   : create RDA files from SL objects
createRDAfiles_fromSLobjects_d29: run_cvsl_varsets_d29
	Rscript code/createRDAfiles_fromSLobjects.R Day29 

## tables_figures   : create tables and figures
tables_figures_d29: createRDAfiles_fromSLobjects_d29
	Rscript code/tables_figures.R Day29
	
######################################################################
## DAY 57 MARKERS ANALYSES
run_cvsl_varsets_d57:
	for number in 1 2 3 4 5 6 7 8 9 10; do \
	Rscript code/run_cvsl_varsets.R $$number Day57; \
	done
    
## createRDAfiles_fromSLobjects   : create RDA files from SL objects
createRDAfiles_fromSLobjects_d57: run_cvsl_varsets_d57
	Rscript code/createRDAfiles_fromSLobjects.R Day57 

## tables_figures   : create tables and figures
tables_figures_d57: createRDAfiles_fromSLobjects_d57
	Rscript code/tables_figures.R Day57
	
######################################################################
## BOTH DAY 29 AND DAY 57 MARKERS ANALYSES
run_cvsl_varsets_Both:
	for number in 1 2 3 4 5 6 7 8 9 10; do \
	Rscript code/run_cvsl_varsets.R $$number Both; \
	done
    
## createRDAfiles_fromSLobjects   : create RDA files from SL objects
createRDAfiles_fromSLobjects_Both: run_cvsl_varsets_Both
	Rscript code/createRDAfiles_fromSLobjects.R Both 

## tables_figures   : create tables and figures
tables_figures_Both: createRDAfiles_fromSLobjects_Both
	Rscript code/tables_figures.R Both
	
######################################################################
## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<
