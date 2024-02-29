#!/usr/bin/env Rscript

## script to compute partial correlation network among phecodes
## Maik Pietzner 06/09/2023
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
require(ppcor)
require(Rfast)

## import phecode to test across proteins
cohort <- args[1]
sex.c  <- args[2]

cat("run partial correlations in", cohort, "among", sex.c ,"\n")

#----------------------------#
##-- import relevant data --##
#----------------------------#

## data set
ukb.dat       <- fread(paste0("input/", cohort))
## subset to sex of interest
if(sex.c != "Both"){
  ukb.dat <- ukb.dat[ sex == sex.c ]
}
## labels
lab.phe       <- fread("input/Labels.phecodes.20230904.txt")
## drop phecodes not covered
lab.phe       <- lab.phe[ id %in% names(ukb.dat) ]

#----------------------------#
##-- compute correlation  --##
#----------------------------#

## compute Pearson correlation
phe.cor           <- Rfast::cora(as.matrix(ukb.dat[, lab.phe$id, with=F]))
# phe.cor       <- cor(as.matrix(ukb.dat[, lab.phe$id[1:10], with=F]))
## add names
rownames(phe.cor) <- colnames(phe.cor) <- lab.phe$id
## convert to table
phe.cor[lower.tri(phe.cor, diag = TRUE)]   <- NA
phe.cor                                    <- na.omit(data.frame(as.table(phe.cor)))

## identify highly correlated pairs
jj                <- subset(phe.cor, Freq > .7)
## delete the one with the larger sample size
jj$drop           <- apply(jj, 1, function(x){
  print(x)
  if(lab.phe$cases.w.primary[which(lab.phe$id == x[1])] >= lab.phe$cases.w.primary[which(lab.phe$id == x[2])]){
    return(x[2])
  }else{
    return(x[1])
  }
})

## drop those from the labels
lab.phe           <- lab.phe[ !(id %in% jj$drop)]

#----------------------------#
##-- compute par. cor.    --##
#----------------------------#

## compute partial correlations
phe.ppcor         <- pcor(ukb.dat[, lab.phe$id, with=F])

## get interesting bits
phe.est            <- phe.ppcor$estimate
rownames(phe.est)  <- colnames(phe.est) <- lab.phe$id
phe.pval           <- phe.ppcor$p.value
rownames(phe.pval) <- colnames(phe.pval) <- lab.phe$id

## convert to data table
phe.est[lower.tri(phe.est, diag = TRUE)]   <- NA
phe.est                                    <- na.omit(data.frame(as.table(phe.est)))
## same for p-values
phe.pval[lower.tri(phe.pval, diag = TRUE)] <- NA
phe.pval                                   <- na.omit(data.frame(as.table(phe.pval)))

## combone both
phe.ppcor          <- merge(phe.est, phe.pval, by=c("Var1", "Var2"), suffixes = c(".estimate", ".pval")) 
## add normal Pearson correlation
phe.ppcor          <- merge(phe.ppcor, phe.cor, by=c("Var1", "Var2"))

## write results to file
fwrite(phe.ppcor, paste("output/partial.correlation", sex.c, cohort, sep="."), sep="\t", row.names=F, na = NA)
