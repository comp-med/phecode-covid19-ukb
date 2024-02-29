#!/usr/bin/env Rscript

## script to run association testing between phecodes and COVID-19 outcomes
## Maik Pietzner 04/09/2023
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

## import type of phecode data set
phe    <- args[1]

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
    tmp  <- ukb.comb[, c("age.covid", "sex.bin", "pop", "non.european", expo, outc, cdat), with=F]
    adj  <- "age.covid + sex.bin + pop"
  }else{
    tmp  <- ukb.comb[ sex == res.surv$sex[x] , c("age.covid", "pop", "non.european", expo, outc, cdat), with=F]
    adj  <- "age.covid + pop"
  }
  
  ## do only if enough observations (>100 events)
  nphe   <- sum(tmp[, ..expo, with=F])
  
  ## run
  if(nphe > 100){
    
    ## where are we at
    cat("run association between", expo, "and", outc, "\n")
    
    ## change adjustment
    adj <- ifelse(outc == "death.covid19", gsub("pop", "non.european", adj), adj)
    
    ## run model: implement cox-prop test and possible step function for time intervals
    ff   <- tryCatch(
      {
        coxph(as.formula(paste0("Surv(", cdat,",", outc, ") ~ ", expo, " + ", adj)), tmp, ties = "breslow")
      }, error=function(e){
        return(NA)
      }, warning=function(w){
        return(NA)
      })
    
    ## proceed only, if model converged
    if(length(ff) > 1){
      
      ## cox-prop test
      ff.p <- cox.zph(ff)
      ## time-varying effect
      ff.s <- survSplit(as.formula(paste0("Surv(", cdat,",", outc, ") ~ ", expo, " + ", adj)), tmp, cut=c(.5,1,1.5,2), episode = "tgroup")
      ## fit the model (need to redefine names of variables)
      ff.s <-  tryCatch(
        {
          coxph(as.formula(paste0("Surv(tstart,", cdat,",", outc, ") ~ ", expo, ":strata(tgroup) + ", adj)), ff.s, ties = "breslow")
        }, error=function(e){
          return(NA)
        })
      
      ## proceed only if model converged
      if(length(ff.s) > 1){
        ## compute summary for storage
        ff   <- summary(ff)
        ff.s <- summary(ff.s)
        ## depend on outcome
        tmp <- data.frame(outcome=outc, phecode=expo, beta=ff$coefficients[1,1], se=ff$coefficients[1,3], pval=ff$coefficients[1,5], nevent=ff$nevent, nall=sum(!is.na(tmp[, outc, with=F])), n.phecode=nphe,
                          p.resid.phecode=ff.p$table[expo, 3], p.resid.overall=ff.p$table["GLOBAL", 3])
        ## add estimates for each strata
        jj  <- grep("strata", rownames(ff.s$coefficients), value=T)
        for(j in 1:length(jj)){
          ## estimate
          tmp[, paste0("beta.", j, ".ti")] <- ff.s$coefficients[jj[j],1]
          ## standard error
          tmp[, paste0("se.", j, ".ti")] <- ff.s$coefficients[jj[j],3]
          ## estimate
          tmp[, paste0("pval.", j, ".ti")] <- ff.s$coefficients[jj[j],5]
        }
        ## return
        return(tmp)
      }
    }
  }
}, mc.cores = 30) 
## combine into one
res.surv <- do.call(plyr::rbind.fill, res.surv)

## write results to file
fwrite(res.surv, paste0("output/Results.", phe), sep = "\t", row.names = F, na = NA)
