#!/bin/bash

echo "Building bookdown $1 report"
Rscript -e "bookdown::clean_book(TRUE)"
Rscript -e "bookdown::render_book(input = 'index_$1.Rmd', config_file = '_bookdown_$1.yml', output_format = bookdown::pdf_document2(keep_tex = TRUE), quiet=TRUE)"
pdflatex -halt-on-error -output-directory _report_$1 _report_$1/covpn_correlates_$1.tex
pdflatex -halt-on-error -output-directory _report_$1 _report_$1/covpn_correlates_$1.tex
pdflatex -halt-on-error -output-directory _report_$1 _report_$1/covpn_correlates_$1.tex
rm _report_$1/*.aux _report_$1/*.log _report_$1/*.lof _report_$1/*.tex _report_$1/*.aux _report_$1/*.toc _report_$1/*.lot