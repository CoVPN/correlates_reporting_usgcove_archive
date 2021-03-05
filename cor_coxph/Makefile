## all           : clean saved results and get results to populate report
all: clean cor_coxph

## clean         : delete saved bootstrap results
clean:
	rm -f output/D29/*
	rm -f output/D57/*

## cor_coxph     : create tables and figures
cor_coxph: dat_seroneg
	Rscript code/cor_coxph.R 57 && Rscript code/cor_coxph.R 29

## dat_seroneg   : create seronegative data set used in analysis
dat_seroneg:
	Rscript code/make_dat_seroneg.R

## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('cor_coxph_report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_coxph
