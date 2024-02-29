#!/bin/sh

## script to perform colocalisation testing for phecodes and COVID-19

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=3:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task=10

#! define how much memory for each node
#SBATCH --mem-per-cpu=2G

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=107,108

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to file>

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
## get paramters
phecode="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F ' ' 'NR == var {print $1}' input/input.phecode.coloc.txt)"
covid="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F ' ' 'NR == var {print $2}' input/input.phecode.coloc.txt)"
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F ' ' 'NR == var {print $3}' input/input.phecode.coloc.txt)"
pos="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F ' ' 'NR == var {print $4}' input/input.phecode.coloc.txt)"


echo "Run coloc of ${phecode} and ${covid} at ${chr}:${pos}"

## set up R environment
source activate Renv

## run the R script
scripts/09_phecode_coloc.R ${phecode} ${covid} ${chr} ${pos}

## deactivate conda env
conda deactivate
