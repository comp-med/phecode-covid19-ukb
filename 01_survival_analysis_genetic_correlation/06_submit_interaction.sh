#!/bin/sh

## script to perform interaction testing for phecode - COVID-19 data

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
#SBATCH --cpus-per-task=32

#! define how much memory for each node
#SBATCH --mem-per-cpu=3G

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-2%10

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to file>

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
## get paramters
cohort="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F ' ' 'NR == var {print $1}' input/interaction.terms.revision)"
inter="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F ' ' 'NR == var {print $2}' input/interaction.terms.revision)"

echo "Run interaction in ${cohort} for ${inter}"

## set up R environment
eval "$(/opt/conda/bin/conda shell.bash hook)"
source activate Renv

## run the R script
scripts/09_interaction_testing.R ${cohort} ${inter}

## deactivate conda env
conda deactivate
