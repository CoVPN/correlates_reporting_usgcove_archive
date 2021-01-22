## data             : clean up and reformat input data, and save to data_clean/
data: code/make_clean_data.R
	Rscript code/make_clean_data.R

## tables           : make summary tables and save data to output/
tables: code/make_table_all.R data
	Rscript code/make_table_all.R

## report           : use tables in output/ to build summary report
report:
	Rscript "rmarkdown::render(input = 'report.Rmd', output_dir = 'output')"

## clean            : delete contents of data_clean/ and output/ directories
clean:
	rm -f data_clean/*.Rdata
	rm -f output/*.Rdata
	rm -f output/*.pdf

## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean data tables
