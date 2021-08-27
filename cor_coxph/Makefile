## all           : clean saved results and get results to populate report
all: clean cor_coxph 

## clean         : delete saved bootstrap results
clean:
	rm -rf output/ $(TRIAL)/*

## cor_coxph     : create tables and figures
cor_coxph:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/cor_coxph.R PrimaryD57
	Rscript code/cor_coxph.R PrimaryD29
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_pooled_real))
	Rscript code/cor_coxph.R PrimaryD29
	#Rscript code/cor_coxph.R PrimaryD29a
endif


## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('cor_coxph_report.Rmd', output_file='cor_coxph_report_$(TRIAL).pdf')"

# turn down number of bootstraps and permutations in CoR for faster builds
before_script:
	echo "B <- 100" >> code/params.R
	echo "numPerm <- 50" >> code/params.R

test: before_script cor_coxph report


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_coxph
