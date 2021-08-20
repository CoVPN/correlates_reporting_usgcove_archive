## all           : clean saved results and get results to populate report
all: clean run_cvsl_varsets createRDAfiles_fromSLobjects tables_figures 

## clean         : delete all SL results
clean:
	rm -f output/*
	rm -f figs/*


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
	Rscript code/run_cvsl_varsets.R 12 && \
	Rscript code/run_cvsl_varsets.R 13 && \
	Rscript code/run_cvsl_varsets.R 14 && \
	Rscript code/run_cvsl_varsets.R 15 && \
	Rscript code/run_cvsl_varsets.R 16 && \
	Rscript code/run_cvsl_varsets.R 17 && \
	Rscript code/run_cvsl_varsets.R 18 && \
	Rscript code/run_cvsl_varsets.R 19 && \
	Rscript code/run_cvsl_varsets.R 20 && \
	Rscript code/run_cvsl_varsets.R 21 && \
	Rscript code/run_cvsl_varsets.R 22 && \
	Rscript code/run_cvsl_varsets.R 23 && \
	Rscript code/run_cvsl_varsets.R 24 && \
	Rscript code/run_cvsl_varsets.R 25 && \
	Rscript code/run_cvsl_varsets.R 26 && \
	Rscript code/run_cvsl_varsets.R 27 && \
	Rscript code/run_cvsl_varsets.R 28 
	
    
## createRDAfiles_fromSLobjects   : create RDA files from SL objects
createRDAfiles_fromSLobjects: run_cvsl_varsets
	Rscript code/createRDAfiles_fromSLobjects.R  

## tables_figures   : create tables and figures
tables_figures: createRDAfiles_fromSLobjects
	Rscript code/tables_figures.R 
	
## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<
