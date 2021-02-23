## all            : clear current results and make figures
all: clean figs cor_boxplots cor_wrcdf

## clean_data     : saves data sets for plotting in data_clean
clean_data: code/data_prep.R
	Rscript code/data_prep.R

## figs           : make figures
figs: fig_code clean_data
	Rscript code/cor_graphics_violin_scatter.R

fig_generation: fig_code data_prep common
    Rscript code/cor_graphics_violin_scatter.R


## boxplots by case and non-case
cor_boxplots: code/cor_data_preprocess.R \
  plot_data params
	Rscript code/CoR_assay_graphics.R

## weighted rcdf: 
cor_wrcdf: code/cor_data_preprocess.R \
  plot_data params
	Rscript code/CoR_wrcdf_with_VE_lines.R

## plot_data        : relabel and clean data sets for plotting
# creates clean_data/long_twophase_data.rds and clean_data/twophase_data.rds
plot_data: code/cor_data_preprocess.R params
	Rscript code/cor_data_preprocess.R

# parameters that may vary across runs
params: code/params.R

## clean          : clean up data_clean and output folder
clean: 
	rm -f data_clean/* && rm -f output/*

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: fig_code figs clean all clean_data

