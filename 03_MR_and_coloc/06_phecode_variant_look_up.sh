#!/bin/sh

## script to extract SNPs

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
#SBATCH --cpus-per-task 2

#! define how much memory for each node
#SBATCH --mem-per-cpu=3G

#! run as job array
#SBATCH --array=1-1445%30

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path to file>

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
## get paramters
file="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F ' ' 'NR == var {print $1}' phecode.files)"

## what is done
echo "look up ${file}"

## grep what is needed
zgrep -wF -f phecode.variants.lookup.20230913.txt <path to file>/${file} > ../tmpdir/phecode.${file}.variant.lookup.results
