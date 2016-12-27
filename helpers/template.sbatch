#!/bin/bash
#################
# This is a template batch script for submission to the Sherlock Cluster at
# Stanford University: http://sherlock.stanford.edu
# It is used by the `mksbatch.sh` script which generates and submits multiple
# batch files.
#################
#
#################
# email all events
#SBATCH --mail-type=ALL
#################
# set a job name
#SBATCH --job-name=XXX_YYY
################
# a file for job output, you can check job progress
#SBATCH --output=XXX_YYY.out
#################
# a file for errors from the job
#SBATCH --error=XXX_YYY.err
#################
# time you think you need; default is one hour
# in minutes in this case, hh:mm:ss
#SBATCH --time=HOURS:00:00
#################
# quality of service; think of it as job priority
#SBATCH -p iric
#################
# number of nodes you are requesting
#SBATCH --nodes=1
#################
# memory per node; default is 4000 MB per CPU
# --mem=4000
#################
# tasks to run per node; a "task" is usually mapped to a MPI processes.
# for local parallelism (OpenMP or threads), use "--ntasks-per-node=1
# --cpus-per-task=16" instead
#SBATCH --ntasks-per-node=1 --cpus-per-task=16
#################

# load modules
module load hdf5
module load gsl/2.1

# run and tee output (redundant)
./RUNFILE | tee out_XXX_YYY
