sbatch 

sbatch tables_figures.sh

results/cvsl_risk_placebo_cvaucs.rda: createRDAfiles_fromSLobjects.sh results/cvsl_riskscore_cvaucs.rds
	sbatch createRDAfiles_fromSLobjects.sh

results/cvsl_riskscore_cvaucs.rds: submit_cluster_job.sh 
	sbatch submit_cluster_job.sh
