## all           : clean saved results and get results to populate report
all: clean cor_nonlinear

## clean         : delete saved bootstrap results
clean:
	rm -f output/*

## cor_nonlinear : create output figures
cor_nonlinear: 
	cd code; Rscript cor_nonlinear.R

## report        : create a standalone report
report:  
	Rscript -e "rmarkdown::render('cor_nonlinear_report.Rmd')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_nonlinear