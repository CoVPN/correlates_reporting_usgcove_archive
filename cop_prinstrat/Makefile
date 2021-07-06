## all           : clean saved results
all: clean cop_prinstrat_covariateAdjusted_condVE 

## clean         : delete saved results
clean:
	rm -rf output/*

## cop_prinstrat_covariateAdjusted_condVE     : run analysis and output results with number of levels of BIP = 4 for example.
cop_prinstrat_covariateAdjusted_condVE:
	Rscript code/cop_prinstrat_covariateAdjusted_condVE.R 4 
ifeq ($(TRIAL), moderna_mock)
	echo "The example code is used to perform analysis on Janssen mock data only"
else ifeq ($(TRIAL), moderna_real)
	echo "The example code is used to perform analysis on Janssen mock data only"
endif

test: clean cop_prinstrat_covariateAdjusted_condVE


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cop_prinstrat_condVE
