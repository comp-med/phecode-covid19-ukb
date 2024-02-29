#!/bin/sh

#BSUB -e <path to file>
#BSUB -o /<path to file>

## script to run the first of two steps for REGENIE
## Maik Pietzner                                  19/10/2021

## force LSF to allocate all cores on one node; implies
## that no job can be submitted if it tries to allocate more cores than a single
## node can offer
#BSUB -R "span[hosts=1]" 

## Request 10GB of memory/core
#BSUB -R "rusage[mem=10GB]"

## how many threads
#BSUB -n 32

## name of the job (do as an array)
#BSUB -J regenie_step1[2-13]%3

## to which queue to be send (type bqueues to see all)
#BSUB -q normal

## get the chromosome
echo "Job ID: $LSB_JOBINDEX"
export batch=${LSB_JOBINDEX}

## export location of genotype files to be used
export dir=<path to file>

## change to relevant directory
cd <path to file>

## run REGENIE (needs a remove flag for samples to be exclude)
<path to file>/regenie_v2.2.4.gz_x86_64_Centos7_mkl \
--step 1 \
--bed ${dir}/ukb_cal_allChrs \
--extract ${dir}/qc_pass.snplist \
--keep ${dir}/qc_pass_whiteEuropean.id \
--phenoFile input/phenotype.${batch}.txt \
--covarFile input/covariates.txt \
--threads 30 \
--bt \
--bsize 1000 \
--lowmem \
--lowmem-prefix tmpdir/regenie_tmp_preds_${batch} \
--out input/ukb_step1_${batch}
