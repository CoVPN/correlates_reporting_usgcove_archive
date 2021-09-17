## all              : tables_output
all: tables_output

## tables_output    : make summary tables and save data to output/
tables_output: code/make_table_all.R 
	Rscript code/make_table_all.R

## clean            : delete contents of data_clean/ and output/ directories
clean:
	rm -f data_clean/*.Rdata
	rm -f output/*.Rdata
	rm -f output/*.pdf

## style            : enforce R code styling conventions automatically
style:
	Rscript -e "styler::style_file('report.Rmd')"
	Rscript -e "styler::style_dir('code')"

## report        : create a standalone report
report:

## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean style tables_output
