## all              
all: clean cor_nonlinear1


## delete saved bootstrap results
clean:
	rm -f output/*


## cor_nonlinear         
cor_nonlinear: cor_nonlinear1 report

cor_nonlinear1: 
	cd code; Rscript cor_nonlinear.R

## standalone report
report:  
	Rscript -e "rmarkdown::render('cor_nonlinear_report.Rmd')"


## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_nonlinear
