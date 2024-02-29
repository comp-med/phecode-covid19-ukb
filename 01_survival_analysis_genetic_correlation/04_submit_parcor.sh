#!/bin/sh

## script to create partial correlation networks for phecodes

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=48:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task=2

#! define how much memory for each node
#SBATCH --mem-per-cpu=20G

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-3

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to file>

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
## get paramters
cohort="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F ' ' 'NR == var {print $1}' input/partial.correlation.cohorts.txt)"
sex="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F ' ' 'NR == var {print $2}' input/partial.correlation.cohorts.txt)"


echo "Run partial correlation in ${cohort} for ${sex}"

## set up R environment
source activate Renv

## run the R script
scripts/07_compute_partial_corr.R ${cohort} ${sex}

## deactivate conda env
conda deactivate
