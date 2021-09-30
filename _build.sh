#!/bin/bash

echo "Building bookdown $1 report"

if [[ "$TRIAL" == "moderna_real" ]] 
then
	if grep -q "predate}" index_cor.Rmd
	then
		echo "Confidentiality statement is already included in index."
	else
		awk '/titling}/{print;print "   - \\predate{\\begin{center}}";print"   - \\postdate{\\vspace{2in} \\\\ \\Large HHS and Moderna Confidential Information – For Official Use Only – Not to be Disseminated INFORMATION NOT RELEASABLE TO THE PUBLIC UNLESS AUTHORIZED BY LAW. This information has not been publicly disclosed and may be a privileged, confidential, deliberative, and/or predecisional communication. It is for internal government use only and must not be disseminated, distributed, or copied to persons not authorized to receive the information. Unauthorized disclosure may result in prosecution to the full extent of the law.\\end{center}}";print;next}1' index_$1.Rmd > tmp
		mv index_$1.Rmd index_$1_tmp.Rmd
		mv tmp index_$1.Rmd
	fi
fi

Rscript -e "bookdown::clean_book(TRUE)"
Rscript -e "bookdown::render_book(input = 'index_$1.Rmd', output_file = 'covpn_correlates_$1_$TRIAL.pdf', config_file = '_bookdown_$1.yml', output_format = bookdown::pdf_document2(keep_tex = TRUE), quiet=TRUE)"
pdflatex -halt-on-error -output-directory _report_$1 _report_$1/covpn_correlates_$1_$TRIAL.tex
pdflatex -halt-on-error -output-directory _report_$1 _report_$1/covpn_correlates_$1_$TRIAL.tex
pdflatex -halt-on-error -output-directory _report_$1 _report_$1/covpn_correlates_$1_$TRIAL.tex
rm -f _report_$1/*.aux _report_$1/*.log _report_$1/*.lof _report_$1/*.tex _report_$1/*.aux _report_$1/*.toc _report_$1/*.lot _report_$1/*.out

if [[ "$TRIAL" == "moderna_real" ]] 
then
	if grep -q "predate}" index_$1.Rmd
	then
		echo "Confidentiality statement is already included in index."
	else
		mv index_$1_tmp.Rmd index_$1.Rmd	
	fi
fi
