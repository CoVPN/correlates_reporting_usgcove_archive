SuperLearner Code
=================

There are 8 important files in this folder:

1.  run\_cvsl\_riskscore.sh (runs R file run\_cvsl\_riskscore.R): This runs CV-Superlearner on the placebo arm and stores SL objects with CV-AUCs in results folder. submit\_cluster\_job.sh is used to run this script by submitting jobs to the cluster.

2.  createRDAfiles\_fromSLobjects.sh (runs createRDAfiles\_fromSLobjects.R): This reads in the SL objects and summarizes the CV-AUCs in an easy-to-read R object.  

3.  tables\_figures.sh (runs tables\_figures.R): This creates SL tables and figures and stores them in output folder. It also stores CV-predictions on the placebo arm. 

4.  constructSL\_predict\_on\_vaccine.sh (runs constructSL\_predict\_on\_vaccine.R): This runs Superlearner, runs predictions on vaccine arm and stores them in output folder.

5. get\_SLweights\_Modelpredictors.sh (runs get\_SLweights\_Modelpredictors.R): This pulls in SL weights and retrieves coefficients of predictors in top-weighted learners. It stores results in output folder.

4.  sl\_screens.R: Contains all screening functions used by the Superlearner.

5.  utils.R: Contains all auxiliary functions used by Superlearner and required in the processing of results.

6.  make\_forest\_plot.R: Contains functions to create forest plots.

The Pipeline of Using the Code
==============================

sbatch run\_cvsl\_riskscore.sh

sbatch createRDAfiles\_fromSLobjects.sh

sbatch tables\_figures.sh

sbatch constructSL\_predict\_on\_vaccine.sh

sbatch get\_SLweights\_Modelpredictors.sh


