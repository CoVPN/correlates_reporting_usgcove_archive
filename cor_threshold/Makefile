## all              : run threshold analysis and get plots
all: clean params data threshold plot

allbutclean: params data threshold plot

## threshold        : run threshold analysis 
threshold: code/Run_Threshold_analysis.R \
  params  
	Rscript code/Run_Threshold_analysis.R

## data             : prepare data for analysis
data: code/clean_data.R params
	Rscript code/clean_data.R
	
## plot             : plot results
plot: code/plotting.R params
	Rscript code/plotting.R

# parameters that may vary across runs
params: code/params.R 

## clean            : delete contents of data_clean/ and figs/
clean:
	rm -f data_clean/*.csv
	rm -f data_clean/Thresholds_by_marker/*.csv
	rm -f figs/simultaneous_CI/*.pdf
	rm -f figs/pointwise_CI/*.pdf
	rm -f figs/pointwise_CI/*.rds
	rm -f figs/simultaneous_CI/*.rds

## style            : enforce R code styling conventions automatically
style:
	Rscript -e "styler::style_file('report.Rmd')"
	Rscript -e "styler::style_dir('code')"

## report        : create a standalone report
report:

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean plot params threshold data
