## all           : clean saved results and get results to populate report
all: clean cor_coxph 

## clean         : delete saved bootstrap results
clean:
	rm -rf output/$(TRIAL)/*

## cor_coxph     : create tables and figures
cor_coxph:
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/cor_coxph.R D29
	Rscript code/cor_coxph.R D57
	Rscript code/cor_graphical_time_to_event.R D29D57
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_pooled_real janssen_na_mock janssen_na_real janssen_la_mock janssen_la_real janssen_sa_mock janssen_sa_real))
	Rscript code/cor_coxph.R D29
	Rscript code/cor_coxph.R D29start1
endif


## report        : create a standalone report
report:
	Rscript -e "rmarkdown::render('cor_coxph_report.Rmd', output_file='cor_coxph_report_$(TRIAL).pdf')"


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_coxph
