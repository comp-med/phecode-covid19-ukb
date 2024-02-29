#!/usr/bin/env Rscript

## script to run interaction association testing between phecodes and COVID-19 outcomes
## Maik Pietzner 18/09/2023
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
require(rms)
require(survival)

## import type of phecode data set (w/ or w/o primary care data)
phe    <- args[1]
int    <- args[2]

cat("run regression with", phe, "\n")

#----------------------------#
##-- import relevant data --##
#----------------------------#

## covariates
ukb.comb <- fread(paste0("input/", phe))
## labels
lab.phe  <- fread("input/Labels.phecodes.20230904.txt")
## drop phecodes not covered
lab.phe  <- lab.phe[ id %in% names(ukb.comb) ]
## keep only variables not included in the interaction
lab.phe  <- lab.phe[ !(id %in% int) ]

#----------------------------#
##--     run testing      --##
#----------------------------#

## do in parallel
registerDoMC(30)

## create grid
res.surv <- merge(as.data.frame(lab.phe[, c("id", "sex")]), data.frame(covid=c("hosp.covid19", "severe_resp", "death.covid19", "long.covid.status"), 
                                                                       date=c("fol.hosp", "fol.covid", "fol.death", "fol.long.covid")))

## run testing
res.surv <- mclapply(1:nrow(res.surv), function(x){
  
  ## get exposure
  expo <- res.surv$id[x]
  ## get outcome
  outc <- res.surv$covid[x]
  ## get date column
  cdat <- res.surv$date[x]
  
  ## date set needed (depends on disease)
  if(res.surv$sex[x] == "Both"){
    tmp  <- ukb.comb[, c("age.covid", "sex", "pop", "non.european", int, expo, outc, cdat), with=F]
    adj  <- paste0("*", paste(unique(c(int, "sex", "age.covid", "non.european")), collapse = " + "))
  }else{
    tmp  <- ukb.comb[ sex == res.surv$sex[x] , c("age.covid", "sex", "pop", "non.european", int, expo, outc, cdat), with=F]
    adj  <- paste0("*", paste(unique(c(int, "age.covid", "non.european")), collapse = " + "))
  }
  
  ## drop age term if needed
  if(int == "age.cat"){
    adj <- gsub(" + age.covid", "", adj, fixed = T)
  }
  
  ## do only if enough observations (>100 events)
  nphe   <- sum(tmp[, ..expo, with=F])
  tphe   <- table(unlist(tmp[, ..expo, with=F]), unlist(tmp[, ..int, with=F]))
  
  ## run
  if(sum(tphe > 50) > 3 & ncol(tphe) > 1){
    
    ## where are we at
    cat("run association between", expo, "and", outc, "with", int, "\n")
    
    ## run model: implement cox-prop test and possible step function for time intervals
    ff   <- tryCatch(
      {
        coxph(as.formula(paste0("Surv(", cdat,",", outc, ") ~ ", expo, adj)), tmp, ties = "breslow")
      }, error=function(e){
        return(NA)
      }, warning=function(w){
        return(NA)
      })
    
    ## proceed only, if model converged
    if(length(ff) > 1){
      
      ## association in each strata
      str  <- unique(unlist(tmp[, ..int]))
      
      ## first level
      ff.1 <- summary(coxph(as.formula(paste0("Surv(", cdat,",", outc, ") ~ ", expo, gsub(paste0("\\*", int), "", adj))), tmp[ eval(as.name(int)) == str[1],], ties = "breslow"))
      ## second level
      ff.2 <- summary(coxph(as.formula(paste0("Surv(", cdat,",", outc, ") ~ ", expo, gsub(paste0("\\*", int), "", adj))), tmp[ eval(as.name(int)) == str[2],], ties = "breslow"))
      
      ## cox-prop test
      ff.p <- cox.zph(ff)
      ## compute summary for storage
      ff   <- summary(ff)
      ## how many entries
      nn   <- nrow(ff$coefficients)
      
      ## depend on outcome
      tmp  <- data.frame(outcome=outc, phecode=expo, 
                         ## main model
                         beta.phecode=ff$coefficients[1,1], se.phecode=ff$coefficients[1,3], pval.phecode=ff$coefficients[1,5], 
                         beta.term=ff$coefficients[2,1], se.term=ff$coefficients[2,3], pval.term=ff$coefficients[2,5],
                         beta.inter=ff$coefficients[nn,1], se.inter=ff$coefficients[nn,3], pval.inter=ff$coefficients[nn,5], 
                         ## stratified models
                         beta.1=ff.1$coefficients[1,1], se.1=ff.1$coefficients[1,3], pval.1=ff.1$coefficients[1,5],
                         beta.2=ff.2$coefficients[1,1], se.2=ff.2$coefficients[1,3], pval.2=ff.2$coefficients[1,5],
                         ## generic information
                         nevent=ff$nevent, nall=sum(!is.na(tmp[, outc, with=F])), n.phecode=nphe,
                         p.resid.phecode=ff.p$table[expo, 3], p.resid.overall=ff.p$table["GLOBAL", 3])
      ## return
      return(tmp)
    }
  }
}, mc.cores = 30) 
## combine into one
res.surv <- do.call(plyr::rbind.fill, res.surv)

## write results to file
fwrite(res.surv, paste0("output/Results.", int, ".", phe), sep = "\t", row.names = F, na = NA)
