#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here("..", "..", ".."))
source(here::here("..", "..", "..", "_common.R"))
#-----------------------------------------------

# this is code from old specs that may be re-implemented for analyses
# other than moderna
if(FALSE){
	# For reproducibility of CV-Superlearner results, use seed “20210217” 
	# and generate 10 seeds by round(runif(10, 1000, 10000)) to run 10 random CV-SL runs.
	set.seed(20210217)
	seeds <- round(runif(10, 1000, 10000))
}else{
	# if only one seed, use this one
	seeds <- 20210216
}

saveRDS(seeds, file = here::here("output", "seeds.rds"))