common_code: ../_common.R
data_prep_code: data_clean/data_prep.R
fig_code: code/cor_graphics_violin_scatter.R
fig_generation: fig_code data_prep common
    Rscript code/cor_graphics_violin_scatter.R



all: cor_boxplots cor_wrcdf

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

## clean            : delete contents of data_clean/ and figs/
clean:
	rm -f data_clean/*.rds
	rm -f figs/*.png
	rm -f Rplots.pdf
