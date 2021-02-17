## all           : clean saved results and get results to populate report
all: clean cor_coxph

## clean         : delete saved bootstrap results
clean:
	rm -f output/*

## cor_coxph     : create tables and figures
cor_coxph: 
	cd code; Rscript cor_coxph.R

## report        : create a standalone report
report:  
	Rscript -e "rmarkdown::render('cor_coxph_report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_coxph