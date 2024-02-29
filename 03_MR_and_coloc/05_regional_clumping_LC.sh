#!/bin/sh

## script to runs regional clumping for LongCOVID statistics

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=2:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task 10

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
  
  ## get the output to be collated
echo "Job ID: $SLURM_ARRAY_TASK_ID"

## add some numbers
export job=$(( $SLURM_ARRAY_TASK_ID ))

## get the output to be collated
echo "Job ID: $job"

## extract the batches to be run
# out="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' input/Outcomes.txt)"
out="st_cas_vs_bd_con"

echo "Collating results for ${out}"

## change to output directory
cd 01_stats

#------------------------------------#
##-- obtain regional lead signals --##
#------------------------------------#

## get all signals of interest (treat extended MHC region as one)
cat ${out}.hg19 | awk '{if(NR == 1) print $0; else if (NR > 1 && ($9+0) < 5e-8) print $0}' - | sed '1d'  | awk '{print $2,$3,$3}' OFS="\t"  | sort -k1,1 -k2,2n | awk '{if($1+0.0==6+0 && $2+0 >= 25500000+500000 && $2+0 <= 34000000-500000) {print $1,25500000,34000000} else if($1+0.0==6+0 && $2+0 >= 25500000 && $2+0 <= 25500000+500000) {print $1,$2-500000,34000000} else if($1+0.0==6+0 && $2+0 >= 34000000-500000 && $2+0 <= 34000000) {print $1,25500000,$2+500000}  else if($2-500000 >= 0) {print $1,$2-500000,$2+500000} else {print $1,0,$2+500000} }' OFS="\t" | /sc-projects/sc-proj-computational-medicine/programs/bedtools merge -i stdin > ../02_regional_signals/regional_results/${out}_regions.txt

## now add the strongest signal for each region
while read chr start end; do cat ${out}.hg19 | awk '{if(NR==1) print $0; else if(NR>1 && ($9+0.0) < 5e-8) print $0}' OFS="\t" - | sed '1d' | awk -v chr=${chr} -v start=${start} -v end=${end} '{if(($2+0)==(chr+0) && ($3+0) >= (start+0) && ($3+0) <= (end+0)) print}' | awk -v start=${start} -v end=${end} -v max=0 '{if(sqrt(($4/$5)*($4/$5))>=max){want=$0; max=sqrt(($4/$5)*($4/$5))}} END{print want,start,end}' OFS="\t" ; done < ../02_regional_signals/regional_results/${out}_regions.txt > ../02_regional_signals/regional_results/${out}_regional_sentinels.txt
