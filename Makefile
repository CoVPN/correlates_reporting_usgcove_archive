##
immuno:
	$(MAKE) -k -C immuno_tabular clean report
	$(MAKE) -k -C immuno_graphical clean report

##
cor:
	$(MAKE) -k -C cor_graphical clean report
	$(MAKE) -k -C cor_coxph clean report
	$(MAKE) -k -C cor_threshold clean report
	$(MAKE) -k -C cor_surrogates clean report
	$(MAKE) -k -C cor_nonpar clean report

##
cop:
	$(MAKE) -k -C cop_prinstrat clean report
	$(MAKE) -k -C cop_controlled clean report
	$(MAKE) -k -C cop_stochastic clean report
	$(MAKE) -k -C cop_mediation clean report

##
style:
	Rscript -e "styler::style_dir(filetype = 'rmd')"

##
book:
	sh ./_build.sh

## type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: style book
