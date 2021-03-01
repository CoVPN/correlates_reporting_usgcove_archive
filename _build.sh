#!/bin/bash

echo "Building bookdown $1 report"
Rscript -e "bookdown::clean_book(TRUE)"
Rscript -e "bookdown::render_book(input = 'index.Rmd', config_file = '_bookdown_$1.yml', output_format = 'bookdown::pdf_document2', quiet=TRUE)"
