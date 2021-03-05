## immuno_analysis        : builds immunogenicity exploratory analyses
immuno_analysis: data_processed
	$(MAKE) -k -C immuno_tabular all
	$(MAKE) -k -C immuno_graphical all

## immuno_report          : builds the CoVPN immunogenicity report
immuno_report: immuno_analysis
	bash ./_build.sh immuno

## cor_analysis           : builds Correlates of Risk analyses
cor_analysis: data_processed
# !!! TMP !!!	$(MAKE) -k -C cor_graphical all
	$(MAKE) -k -C cor_coxph all
# these will be added to cor report eventually
#	$(MAKE) -k -C cor_threshold all
#	$(MAKE) -k -C cor_surrogates all
#	$(MAKE) -k -C cor_nonpar all

## cor_report             : builds the CoVPN correlates of risk report
cor_report: cor_analysis
	bash ./_build.sh cor

## cop_analysis           : builds Correlates of Protection analyses
cop_analysis: data_processed
	$(MAKE) -k -C cop_prinstrat all
	$(MAKE) -k -C cop_controlled all
	$(MAKE) -k -C cop_stochastic all
	$(MAKE) -k -C cop_mediation all

## cop_report             : builds the CoVPN correlates of protection report
cop_report: cop_analysis
	bash ./_build.sh cop

## data_processed         : create processed data from raw data
# make_dat_proc.R needs to be executed in data_clean folder to correctly
# active renv
data_processed:
	Rscript data_clean/make_dat_proc.R

## style                  : re-styles the codebase for consistent formatting
style:
	Rscript -e "styler::style_dir(filetype = 'rmd')"

## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: style help immuno_analysis \
  immuno_report cor_report cor_analysis data_processed
