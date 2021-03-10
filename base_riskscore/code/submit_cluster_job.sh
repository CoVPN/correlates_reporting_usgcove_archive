#!/bin/bash
ml fhR/4.0.2-foss-2019b
ml jbigkit
sbatch -c10 --mem 11G --time=7-0 --constraint=gizmok code/run_cvsl_riskscore.sh
