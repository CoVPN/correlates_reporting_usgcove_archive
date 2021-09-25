#!/bin/bash

echo $TRIAL
#Rscript -e "bookdown::clean_book(TRUE)"
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_$1_$TRIAL.pdf', config_file = '_bookdown_$1.yml', output_format = bookdown::pdf_document2(keep_tex = TRUE), quiet=TRUE)"
rm -f _report_cor/*.aux _report_cor/*.log _report_cor/*.lof _report_cor/*.tex _report_cor/*.aux _report_cor/*.toc _report_cor/*.lot _report_cor/*.out
