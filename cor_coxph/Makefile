## all              
all: clean cor_coxph1


## delete saved bootstrap results
clean:
	rm -f output/*


## cor_coxph         
cor_coxph: cor_coxph1 report

cor_coxph1: 
	cd code; Rscript CoR.R

## standalone report
report:  
	Rscript -e "rmarkdown::render('cor_coxph_report.Rmd')"


## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_coxph




