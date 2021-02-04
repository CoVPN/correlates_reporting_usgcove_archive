## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

## immuno_analysis        : builds immunogenicity exploratory analyses
immuno_analysis:
	$(MAKE) -k -C immuno_tabular all
	$(MAKE) -k -C immuno_graphical all

## immuno_report          : builds the CoVPN immunogenicity report
immuno_report:
	bash ./_build.sh immuno

## cor_analysis           : builds Correlates of Risk analyses
cor_analysis:
	$(MAKE) -k -C cor_graphical all
	$(MAKE) -k -C cor_coxph all
	$(MAKE) -k -C cor_threshold all
	$(MAKE) -k -C cor_surrogates all
	$(MAKE) -k -C cor_nonpar all

## cor_report             : builds the CoVPN correlates of risk report
cor_report:
	bash ./_build.sh cor

## cop_analysis           : builds Correlates of Protection analyses
cop_analysis:
	$(MAKE) -k -C cop_prinstrat all
	$(MAKE) -k -C cop_controlled all
	$(MAKE) -k -C cop_stochastic all
	$(MAKE) -k -C cop_mediation all

## cop_report             : builds the CoVPN correlates of protection report
cop_report:
	bash ./_build.sh cop

## style                  : re-styles the codebase for consistent formatting
style:
	Rscript -e "styler::style_dir(filetype = 'rmd')"

.PHONY: style immuno_analysis cor_analysis cop_analysis
