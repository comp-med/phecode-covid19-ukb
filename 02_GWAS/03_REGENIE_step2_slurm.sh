#!/bin/sh

## script to runs econd step of REGENIE

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
#SBATCH --cpus-per-task 120

#! no other jobs can be run on the same node
#SBATCH --exclusive

#! run as job array
#SBATCH --array=1-6%5

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## export location of genotype files to be used
export dir=<path to file>
## directory of variant inclusion files
export var=<path to file>

## get the chromosome
echo "Job ID: $SLURM_ARRAY_TASK_ID"

## extract the batches to be run
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' input/Job.submission.REGENIE.step2.phecodes.rerun.txt)"
batch="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' input/Job.submission.REGENIE.step2.phecodes.rerun.txt)"

##change to X for 23
if [[ $chr -eq 23 ]]; then

  export chr="X"
  
  cd <path to file>
  
  ## create variant inclusion list for the respective chromosome
  cat ${var}/output/ukb_imp_chr${chr}_qced.txt | awk -v chr=${chr} '{if(NR != 1 && $3 == chr) print $2}' - > tmpdir/tmp.ex.chr${chr}.${batch}.list
  
  ## run REGENIE (needs a remove flag for samples to exclude)
  <path to file>/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
  --step 2 \
  --bgen ${dir}/ukb_imp_chr${chr}_v3.bgen \
  --ref-first \
  --extract tmpdir/tmp.ex.chr${chr}.${batch}.list \
  --sample  ${dir}/ukb_imp_chrX_v3.sample \
  --keep ${var}/input/UKBB.samples.51157.inclusion.white.20210927.sample \
  --phenoFile input/phenotype.${batch}.txt \
  --covarFile input/covariates.txt \
  --threads 120 \
  --bt \
  --spa 0.01 \
  --pred input/ukb_step1_${batch}_pred_new.list \
  --bsize 800 \
  --out output/ukb_step2_${batch}_chr${chr}
  
  ## remove variant inclusion list
  rm tmpdir/tmp.ex.chr${chr}.${batch}.list

else
  
  cd <path to file>

  ## create variant inclusion list for the respective chromosome
  cat ${var}/output/ukb_imp_chr${chr}_qced.txt | awk -v chr=${chr} '{if(NR != 1 && $3 == chr) print $2}' - > tmpdir/tmp.ex.chr${chr}.${batch}.list
  
  ## run REGENIE (needs a remove flag for samples to exclude)
  <path to file>/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
  --step 2 \
  --bgen ${dir}/ukb_imp_chr${chr}_v3.bgen \
  --ref-first \
  --extract tmpdir/tmp.ex.chr${chr}.${batch}.list \
  --sample  ${dir}/ukb_imp_chr1_v3.sample \
  --keep ${var}/input/UKBB.samples.51157.inclusion.white.20210927.sample \
  --phenoFile input/phenotype.${batch}.txt \
  --covarFile input/covariates.txt \
  --threads 120 \
  --bt \
  --spa 0.01 \
  --pred input/ukb_step1_${batch}_pred_new.list \
  --bsize 800 \
  --out output/ukb_step2_${batch}_chr${chr}
  
  ## remove variant inclusion list
  rm tmpdir/tmp.ex.chr${chr}.${batch}.list
fi