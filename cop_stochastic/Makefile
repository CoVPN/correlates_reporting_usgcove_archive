## all              : run stochastic interventional analysis and make figures
all: clean data analysis figs

## data             : prepare data for analysis
data: code/make_analysis_data.R
	Rscript code/make_analysis_data.R

## analysis         : run stochastic interventional risk analysis
analysis: code/run_cop_analysis.R
	Rscript code/run_cop_analysis.R

## figs             : make figures of analysis results
figs: code/make_sve_cop_figs.R analysis
	Rscript code/make_sve_cop_figs.R

## clean            : delete contents of data_clean/ and figs/
clean:
	rm -f data_clean/*.rds
	rm -f output/*.rds
	rm -f figs/*.pdf

## style            : enforce R code styling conventions automatically
style:
	Rscript -e "styler::style_file('report.Rmd')"
	Rscript -e "styler::style_dir('code')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

report:
	Rscript -e "renv::activate(here::here('..')); rmarkdown::render('report.Rmd')"

# define PHONYS
.PHONY: help clean params data
