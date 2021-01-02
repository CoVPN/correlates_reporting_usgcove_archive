## report                    : compiles the final report in output/report.html
report: code/report.Rmd \
  output/summ_stats.RData \
  figs/fig1.png figs/fig2.png
	cd code && Rscript make_report.R

## clean_data/clean_data.txt : two raw data sets merged
clean_data/clean_data.txt: code/make_clean_data.R \
  raw_data/biomarkers.txt raw_data/clinical.txt
	Rscript code/make_clean_data.R

## output/summ_stats.txt     : summary statistics reported in text of report
output/summ_stats.RData: code/make_summary_stats.R \
  code/good_round.R \
  clean_data/clean_data.txt
	Rscript code/make_summary_stats.R

## figs/fig1.png             : scatter plot of biomarker 1 vs. biomarker
figs/fig1.png: code/make_fig1.R clean_data/clean_data.txt
	Rscript code/make_fig1.R

## figs/fig1.png             : paneled scatter plot of biomarkers by age, sex
figs/fig2.png: code/make_fig2.R clean_data/clean_data.txt
	Rscript code/make_fig2.R

## figs                      : populates /figs with all figures
.PHONY: figs
figs: figs/fig1.png figs/fig2.png

## clean                     : removes all contents of output/ and figs/
clean:
	rm clean_data/* figs/* output/*

help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: report clean help
