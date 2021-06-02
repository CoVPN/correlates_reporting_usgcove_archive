## all           : clean saved results and get results to populate report
all: clean cor_nonlinear

## clean         : delete saved bootstrap results
clean:
	rm -rf output/*

## cor_nonlinear : create output figures
cor_nonlinear: 
	Rscript code/cor_nonlinear.R 29
ifeq ($(TRIAL), moderna_mock)
	Rscript code/cor_nonlinear.R 57
endif

## report        : create a standalone report
report:  
	Rscript -e "rmarkdown::render('cor_nonlinear_report.Rmd', output_file='cor_nonlinear_report_$(TRIAL).pdf')"

# turn down number of bootstraps and permutations in CoR for faster builds
before_script:
	echo "B <- 100" >> code/params.R
	echo "numPerm <- 50" >> code/params.R

test: before_script clean cor_nonlinear report


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_nonlinear
