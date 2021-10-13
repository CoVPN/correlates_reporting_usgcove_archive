## all            : clear current results and make figures
all: clean figs

## clean_data     : saves data sets for plotting in data_clean
clean_data: code/cor_data_preprocess.R
	Rscript code/cor_data_preprocess.R D29D57

## figs           : make figures
figs: fig_code clean_data
	Rscript code/cor_graphics_violin_scatter.R D29D57
	Rscript code/CoR_assay_graphics.R D29D57
	Rscript code/CoR_wrcdf_with_VE_lines.R D29D57

fig_code: code/cor_graphics_violin_scatter.R \
  code/CoR_assay_graphics.R \
  code/CoR_wrcdf_with_VE_lines.R \
  code/params.R

## clean          : clean up data_clean and output folder
clean:
	rm -f data_clean/* && rm -rf output/*

## report        : create a standalone report
report:

# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

.PHONY: fig_code figs clean all clean_data
