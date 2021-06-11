## all              : run threshold analysis and get plots
all: results

## data             : prepare data for analysis
data_clean/data_Day57.rds: code/clean_data.R code/params.R
	Rscript code/clean_data.R
data_clean/data_Day57_cat.rds: data_clean/data_Day57.rds
data_clean/data_Day29.rds: data_clean/data_Day57.rds
data_clean/data_Day29_cat.rds: data_clean/data_Day57.rds	

data: data_clean/data_Day57.rds data_clean/data_Day29.rds \
  data_clean/data_Day57_cat.rds data_clean/data_Day29_cat.rds

## results          : get mediation effects for all assays
output/full_result.rds: code/run_analysis.R data
	Rscript $<
output/full_result_cat.rds: output/full_result.rds
## Recover from the removal of $@
	@test -f $@ || rm -f output/full_result.rds
	@test -f $@ || $(MAKE) $(AM_MAKEFLAGS) output/full_result.rds

results: output/full_result.rds output/full_result_cat.rds

## report           : build a report locally
report: report.Rmd output/full_result.rds output/full_result_cat.rds
	Rscript -e "rmarkdown::render('report.Rmd', output_format = 'pdf_document')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help params data report results