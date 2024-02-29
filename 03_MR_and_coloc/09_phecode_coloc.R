#!/usr/bin/env Rscript

## script to run colocalisation between phecodes and COVID-19 outcomes
## Maik Pietzner 20/09/2023
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)
# print(R.Version())

## correct directory
setwd("<path to file>")

## packages needed
require(data.table)
require(doMC)
require(coloc)
require(susieR)
require(biomaRt)

## import type of phecode data set (w/ or w/o primary care data)
phecode <- args[1]
covid   <- args[2]
chr.s   <- as.numeric(args[3])
pos.m   <- as.numeric(args[4])

cat("run coloc between", phecode, "and", covid, "at", chr.s, pos.m, "\n")

#--------------------------------------------#
##--             import data              --##
#--------------------------------------------#

## augment region to query
pos.s    <- max(c(0,pos.m-5e5))
pos.e    <- pos.m+5e5

## --> phecode <-- ##

## read the relevant data
res.phe            <- paste0("zcat <path to file>/gwas_results/", phecode, ".allchr.results.gz | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                         " '{if(($1 == chr && $2 >= low && $2 <= upp) || NR == 1) print $0}'")
res.phe            <- data.table::fread(cmd = res.phe, sep = " ", header = T)
## drop SNPs that have possibly failed in REGENIE
res.phe            <- res.phe[is.na(EXTRA)]
## add MarkerName
res.phe[, MarkerName := paste0("chr", CHROM, ":", GENPOS, "_", pmin(ALLELE0, ALLELE1), "_", pmax(ALLELE0, ALLELE1))]

## --> COVID-19 <-- ##

## import data
res.cov            <- paste0("zcat <path to file>/", covid, "_ALL_exUKBB.tsv.gz | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                             " '{if(($1 == chr && $2 >= low && $2 <= upp) || NR == 1) print $0}'")
res.cov            <- data.table::fread(cmd = res.cov, sep = "\t", header = T)
## rename
names(res.cov)[1]  <- "CHR"
## keep only what is needed
res.cov            <- res.cov[, c("CHR", "POS", "REF", "ALT", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_inv_var_meta_p")]
## rename
names(res.cov)     <- c("CHROM", "GENPOS", "ALLELE0", "ALLELE1", "BETA", "SE", "PVAL")
## add log10p
res.cov[, LOG10P := -log10(PVAL)]
## create MarkerName
res.cov[, MarkerName := paste0("chr", CHROM, ":", GENPOS, "_", pmin(ALLELE0, ALLELE1), "_", pmax(ALLELE0, ALLELE1))]

## --> combine <-- ##

## combine the two
res.comb           <- merge(res.phe, res.cov, by=c("MarkerName", "CHROM", "GENPOS"), suffixes = c(".phecode", ".covid"))
## align effect estimates
res.comb[, BETA.covid := ifelse(ALLELE1.covid == ALLELE1.phecode, BETA.covid, -BETA.covid)]
## drop what is no longer needed
res.comb$ALLELE0.covid <- res.comb$ALLELE1.covid <- NULL

#-----------------------------------------#
##--            LD-matrix              --##
#-----------------------------------------#

## write list of SNPs to be queried to file 
write.table(res.comb$ID, paste("tmpdir/snplist", paste(phecode, covid, sep="."), chr.s, pos.s, pos.e, "lst", sep="."), col.names = F, row.names = F, quote = F)

## snp data to be queried
tmp.z        <- res.comb[, c("ID", "CHROM", "GENPOS", "ALLELE0.phecode", "ALLELE1.phecode")]
names(tmp.z) <- c("rsid", "chromosome", "position", "allele1", "allele2")

## adopt chromosome if needed
if(chr.s < 10){
  print("tada")
  tmp.z$chromosome <- paste0("0", tmp.z$chromosome)
}else if(chr.s == 23){
  tmp.z$chromosome <- "X"
}

## check for input
print(head(tmp.z))

## write to file
write.table(tmp.z, paste("tmpdir/snpz", paste(phecode, covid, sep="."), chr.s, pos.s, pos.e, "z", sep="."), row.names = F, quote = F)

## --> create master file for LDstore2 <-- ##

## assign entries
m.file      <- data.frame(z=paste("tmpdir/snpz", paste(phecode, covid, sep="."), chr.s, pos.s, pos.e, "z", sep="."),
                          bgen=paste("tmpdir/filtered", paste(phecode, covid, sep="."), chr.s, pos.s, pos.e, "bgen", sep="."),
                          bgi=paste("tmpdir/filtered", paste(phecode, covid, sep="."), chr.s, pos.s, pos.e, "bgen.bgi", sep="."),
                          ld=paste("tmpdir/ld", paste(phecode, covid, sep="."), chr.s, pos.s, pos.e, "ld", sep="."),
                          incl="<path to file>/UKBB.samples.51157.inclusion.white.unrelated.20210927.incl",
                          n_samples=ifelse(chr.s != 23, 487409, 486757))

## write to file
write.table(m.file, paste("tmpdir/master", paste(phecode, covid, sep="."), chr.s, pos.s, pos.e, "z", sep="."), sep=";", row.names = F, quote = F)

cat("--------------------------------------\n")
cat("computing LD matrix\n")

## submit the job
system(paste("./scripts/get_LD_matrix_ldstore.sh", chr.s, pos.s, pos.e, paste(phecode, covid, sep=".")))

## read in matrix
ld             <- fread(paste("tmpdir/ld", paste(phecode, covid, sep="."), chr.s, pos.s, pos.e, "ld", sep="."), data.table = F)
## print for checking
print(ld[1:5,1:5])

## create identifier column in results to keep the mapping to the LD matrix
res.comb$snp.id <- 1:nrow(res.comb) 

## rename
row.names(ld) <- colnames(ld) <- res.comb$MarkerName

## delete large file
system(paste("rm tmpdir/ld", paste(phecode, covid, sep="."), chr.s, pos.s, pos.e, "ld", sep="."))

cat("Done\n")
cat("--------------------------------------\n")

#-----------------------------------------#
##--    prepare coloc (using SuSiE)    --##
#-----------------------------------------#

## compute MAF
res.comb[, MAF := ifelse(A1FREQ < .5, A1FREQ, 1 - A1FREQ)]

## --> prepare input <-- ##

## phecode
D1          <- list(beta=res.comb$BETA.phecode, varbeta=res.comb$SE.phecode^2, 
                    type="cc",
                    MAF=res.comb$MAF,
                    snp=res.comb$MarkerName,
                    position=1:nrow(res.comb),
                    ## subset the LD matrix to what is needed
                    LD=as.matrix(ld[res.comb$MarkerName, res.comb$MarkerName]))

## COVID-19
D2          <- list(beta=res.comb$BETA.covid, varbeta=res.comb$SE.covid^2, 
                    type="cc",
                    MAF=res.comb$MAF,
                    snp=res.comb$MarkerName,
                    position=1:nrow(res.comb),
                    ## subset the LD matrix to what is needed
                    LD=as.matrix(ld[res.comb$MarkerName, res.comb$MarkerName]))

## --> fine-mapping <-- ##

## Phecode
set.seed(42)
susie.phecode <- tryCatch(
  {
    runsusie(D1, max_iter = 10000, L=5)
  }, error=function(e){
    return(NA)
  })

## trait
susie.covid <- tryCatch(
  {
    runsusie(D2, max_iter = 10000, L=5, prior_variance=0.2^2, estimate_prior_variance=FALSE)
  }, error=function(e){
    return(NA)
  })

## --> colocalisation <-- ##

## check whether both have outcome
if(length(susie.phecode) > 1 & length(susie.covid) > 1){
  
  ## additional check, whether both traits have at least one credible set to be tested
  if(length(summary(susie.phecode)$cs) > 0 & length(summary(susie.covid)$cs) > 0){
    
    #-----------------------------------------#
    ##-- 	            run coloc            --##
    #-----------------------------------------#
    
    ## run coloc with susie input
    res.coloc        <- coloc.susie(susie.phecode, susie.covid, p12 = 5e-6)
    
    #-----------------------------------------#
    ##--   cross-check with fine-mapping   --##
    #-----------------------------------------#
    
    ## add LD with fine-mapped variants 
    res.coloc         <- res.coloc$summary
    
    print(res.coloc)
    
    ## compute ld between selected lead hits
    res.coloc$ld.top  <- apply(res.coloc[, c("hit1", "hit2"), drop=F], 1, function(k) ld[k[1], k[2]]^2)
    
    ## sanity check
    res.coloc$ld.sens <- sapply(res.coloc$hit2, function(x) ld[x, res.comb$MarkerName[which(res.comb$GENPOS == pos.m)]]^2)
    
    ## add effect estimates (based on COVID-19 variant)
    res.coloc         <- merge(res.coloc, res.comb, by.x = "hit2", by.y = "MarkerName")
    
    ## --> to ease import <-- ##
    
    ## import label data
    lab.phe                     <- fread("<path to file>/Labels.phecodes.20230904.txt")
    lab.cov                     <- data.frame(id = c("A2", "B2", "C2", "LC"), phenotype = c("COVID-19 critical illness", "COVID-19 hospitalisation", "SARS-CoV2 infection", "Long COVID"))
    ## add label
    res.coloc$id.phecode        <- phecode
    res.coloc$phenotype.phecode <- lab.phe$phenotype[which(lab.phe$id == phecode)]
    res.coloc$id.covid          <- covid
    res.coloc$phenotype.covid   <- lab.cov$phenotype[which(lab.cov$id == covid)]
    ## add original hit
    res.coloc$covid19.lead      <- res.comb$MarkerName[which(res.comb$GENPOS == pos.m)]
    
    #-----------------------------------------#
    ##-- 	        draw selected            --##
    #-----------------------------------------#
    
    if(sum(res.coloc$PP.H4.abf > .7) > 0 | sum(res.coloc$ld.top > .8) > 0){

      ## import plotting function
      source("scripts/plot_locus_compare.R")
      
      ## actual plot
      png(paste("graphics/susie", phecode, covid, chr.s, pos.s, pos.e, "png", sep="."), width=16, height=8, units="cm", res=200)
      par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, cex.main=.6, font.main=2)
      ## more complex layout for gene assignment
      layout(matrix(c(1,1,1,2,3,4),3,2), heights = c(.43,.37,.2))
      plot.locus.compare(res.comb, res.coloc, ld, a.vars=unique(c(res.coloc$hit1, res.coloc$hit2)))
      dev.off()
    }
    
    #-----------------------------------------#
    ##--          	store results          --##
    #-----------------------------------------#
    
    write.table(res.coloc, paste("output/susie", phecode, covid, chr.s, pos.s, pos.e, "txt", sep="."), sep = "\t", row.names = F)
    
  }else{
    
    cat("No results from fine-mapping\n")
    write.table(data.frame(id.phecode=NA), paste("output/susie", phecode, covid, chr.s, pos.s, pos.e, "txt", sep="."), sep = "\t", row.names = F)
    
  }

  }else{
  
    cat("No results from fine-mapping\n")
    write.table(data.frame(id.phecode=NA), paste("output/susie", phecode, covid, chr.s, pos.s, pos.e, "txt", sep="."), sep = "\t", row.names = F)
    
}

