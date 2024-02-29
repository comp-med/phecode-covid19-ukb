#!/bin/sh
## Author Dr. Summaira Yasmeen

## script to run compute genetic correlations
## Requires LDSC: https://github.com/bulik/ldsc
## HPC SLURM
## script to runs Genetic Correlation(rg) between COVID-19 and phecodes 

#SBATCH --job-name Genetic Correlation
#SBATCH --partition=compute
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 120
#SBATCH --array=1-4%4
#SBATCH --mail-type=FAIL
#SBATCH --output=~../log/Genetic Correlation-%j.log



## project in running dir..
export current_proj=.../genetic_correlations
export input=${current_proj}/input
export output=${current_proj}/output

## output 
export tmpdir=${current_proj}/tmpdir
export LDSC=../ldsc

# get covid_19 data 
#Strict case definition (Long COVID after test-verified SARS-CoV-2 infection) vs broad control definition (population control)
curl -o "${input}/st_cas_vs_bd_con"  "https://my.locuszoom.org/gwas/192226/data/?token=09a18cf9138243db9cdf79ff6930fdf8" --compressed
#Broad case definition (Long COVID after any SARS-CoV-2 infection) vs broad control definition
curl -o "${input}/bd_cas_vs_bd_con.gz"  "https://my.locuszoom.org/gwas/826733/data/?token=c7274597af504bf3811de6d742921bc8" --compressed
#Strict case definition vs strict control definition (individuals that had SARS-CoV-2 but did not develop Long COVID)
curl -o "${input}/st_cas_vs_st_con.gz"  "https://my.locuszoom.org/gwas/793752/data/?token=0dc986619af14b6e8a564c580d3220b4" --compressed
#Broad case definition vs strict control definition:
curl -o "${input}/bd_cas_vs_st_con.gz"  "https://my.locuszoom.org/gwas/91854/data/?token=723e672edf13478e817ca44b56c0c068" --compressed

ls ${input}/*.gz | xargs -n1 basename -s .gz > ${input}/covid_19.list.txt
covid_19="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' ${input}/covid_19.list.txt)"

## specify the rg string 
## example rg_string = st_cas_vs_bd_con.gz,phe1_gwas,phe2_gwas,...,pheN_gwas
rg_string="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' ${input}/rg.string)"

echo "Perform Gene_Correlation for ${covid_19} .. Job ID: $SLURM_ARRAY_TASK_ID"

## -->  Munge summary stats <-- ##
${LDSC}/munge_sumstats.py \
	--sumstats ${tmpdir}/${covid_19} \
	--snp SNP \
	--a1 EA \
	--a2 NEA \
	--p Pval \
	--signed-sumstats BETA,0 \
	--out ${tmpdir}/${covid_19} \
	--chunksize 500000 \
	--keep-maf \
	--merge-alleles ${LDSC}/w_hm3.snplist

## --> run gene correlation<-- ##

/ldsc.py \
	--rg  ${rg_string} \
	--ref-ld-chr ${LDSC}/eur_w_ld_chr/ \
	--w-ld-chr ${LDSC}/eur_w_ld_chr/ \
	--out ${output}/${covid_19}_rg


