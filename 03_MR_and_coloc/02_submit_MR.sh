#!/bin/sh

## script to run MR analysis with phecodes as exposure and COVID-19 as outcome

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=72:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task 5

#! define how much memory for each node
#SBATCH --mem-per-cpu=4G

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-49%20

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to file>

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
expo="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' input/mr.input.phecode.covid19.txt)"

echo "Node ID: $SLURM_NODELIST"

echo "Exposure: ${expo}"

## set up R environment
source activate Renv

## run the R script
scripts/03_run_MR_phecode_COVID19.R ${expo}

## deactivate conda env
conda deactivate
