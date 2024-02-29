#!/bin/sh

## script to generate plots for regional diversity

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=12:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task=15

#! define how much memory for each node
#SBATCH --mem-per-cpu=10G

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to file>

## get the exposure (add one to avoid the header)
echo "Job ID: $SLURM_ARRAY_TASK_ID"

## set up R environment
eval "$(/opt/conda/bin/conda shell.bash hook)"
source activate Renv

## run the R script
scripts/11_association_testing_extended.R "UKB.phecodes.COVID19.20240102.txt"

## deactivate conda env
conda deactivate
