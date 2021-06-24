SuperLearner Code
=================

There are 6 files in this folder:

1.  run\_sl\_assays.sh (use to run R file run\_sl\_risk\_score.R)

2.  submit\_sl\_assays\_updated.sh (submit jobs to the cluster, each job
    returns a single superlearner model)

3.  run\_sl\_risk\_score.R (train the superlearner models)

4.  sl\_screens.R (screening functions)

5.  utils.R (axuliary functions)

6.  sl\_prediction\_risk\_scores.R (predict risk scores)

The Pipeline of Using the Code
==============================

Submit a series of jobs by using submit\_sl\_assays\_updated.sh. Here
you can define the number of jobs you want. Each job id represents a
different task and you can go to the run\_sl\_risk\_score.R, there are
some commnets in the last part describing the correspondence of job id.
After submiiting the job and get a number of supearlearner models,
sl\_prediction\_risk\_scores.R can be used to predict the risk scores
