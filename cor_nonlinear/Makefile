## all           : clean saved results and get results to populate report
all: clean cor_nonlinear

## clean         : delete saved bootstrap results
clean:
	rm -rf output/$(TRIAL)/*

## cor_nonlinear : create output figures
cor_nonlinear: 
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/cor_nonlinear.R D57
	Rscript code/cor_nonlinear.R D29
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_pooled_real))
	Rscript code/cor_nonlinear.R D29
	Rscript code/cor_nonlinear.R D29start1
endif

## report        : create a standalone report
report:  
	Rscript -e "rmarkdown::render('cor_nonlinear_report.Rmd', output_file='cor_nonlinear_report_$(TRIAL).pdf')"

test: clean cor_nonlinear report


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_nonlinear
