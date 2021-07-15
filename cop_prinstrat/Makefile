## all           : clean saved results and get results to populate report
all: clean cop_prinstrat 

## clean         : delete saved bootstrap results
clean:
	rm -rf output/*

## cop_prinstrat     : create tables and figures. If TRIAL=janssen pooled mock, run condVE analysis with number of levels of BIP = 4 for example. Number of replications for bootstrap is set to be 500 with number of cores to be detectCores().
cop_prinstrat:
ifeq ($(TRIAL), moderna_mock)
	Rscript code/cop_cove_Run.R 57
	Rscript code/cop_cove_Run.R 29
else ifeq ($(TRIAL), moderna_real)
	Rscript code/cop_cove_Run.R 57
	Rscript code/cop_cove_Run.R 29
else ifeq ($(TRIAL), janssen_pooled_mock)
	Rscript code/cop_ensemble_Run.R 29
	Rscript code/cop_prinstrat_covariateAdjusted_condVE.R 4
else ifeq ($(TRIAL), janssen_pooled_real)
	Rscript code/cop_ensemble_Run.R 29
endif

## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('cop_prinstrat_report.Rmd', output_file='cop_prinstrat_report_$(TRIAL).pdf')"

# turn down number of bootstraps and permutations in cor for faster builds
before_script:
	echo "b <- 100" >> code/params.r
	echo "numperm <- 50" >> code/params.r

test: before_script clean cop_prinstrat report


# type 'make help' to show all make commands
help: makefile
	@sed -n 's/^##//p' $<

# define phonys
.PHONY: help clean cop_prinstrat