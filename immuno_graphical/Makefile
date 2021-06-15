## all              : twophase_plots + demog_plots
all: clean twophase_plots demog_plots

## twophase_plots   : rcdfs, boxplots, scatter plots of assays saved in figs/
twophase_plots: code/descriptive_graphics_two_phase_plots.R \
  params plot_data make_assay_lim
	Rscript code/descriptive_graphics_two_phase_plots.R

## demog_plots      : rcdfs, boxplots by demographics saved in figs/demographics
demog_plots: code/descriptive_graphics_demographic_strata.R \
  params plot_data make_assay_lim
	Rscript code/descriptive_graphics_demographic_strata.R

## plot_data        : relabel and clean data sets for plotting
# creates clean_data/long_twophase_data.rds and clean_data/twophase_data.rds
plot_data: code/descriptive_graphics_data_preprocess.R params
	Rscript code/descriptive_graphics_data_preprocess.R

# parameters that may vary across runs
params: code/params.R

# define the axis limits for the assays
make_assay_lim: params plot_data code/make_assay_lim.R 
	Rscript code/make_assay_lim.R

## clean            : delete contents of data_clean/ and figs/
clean:
	rm -f data_clean/*.rds
	rm -f figs/*.png
	rm -f figs/demographics/*.png
	rm -f Rplots.pdf

## style            : enforce R code styling conventions automatically
style:
	Rscript -e "styler::style_file('report.Rmd')"
	Rscript -e "styler::style_dir('code')"

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean plot_data params make_assay_lim demog_plots twophase_plots
