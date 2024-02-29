#!/usr/bin/env Rscript

## script to run association testing between phecodes and COVID-19 outcomes, using a multi-state model
## with death from any cause as competing risk
## Maik Pietzner 03/01/2024
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
# phe    <- "UKB.phecodes.COVID19.20240102.txt"

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
res.surv <- merge(as.data.frame(lab.phe[, c("id", "sex")]), data.frame(covid=c("hosp.covid19.comprsk", "severe_resp.comprsk", "death.covid19.comprsk", "long.covid.status.comprsk"), 
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
    tmp  <- ukb.comb[, c("f.eid", "age.covid", "sex.bin", "pop", "non.european", expo, outc, cdat), with=F]
    adj  <- "age.covid + sex.bin + pop"
  }else{
    tmp  <- ukb.comb[ sex == res.surv$sex[x] , c("f.eid", "age.covid", "pop", "non.european", expo, outc, cdat), with=F]
    adj  <- "age.covid + pop"
  }
  
  # ## do only if enough observations (>100 events)
  # nphe   <- sum(tmp[, ..expo, with=F])
  
  ## careful, model will not converge, if no death event occurs in people with the disease
  nphe   <- table(unlist(tmp[, ..outc, with=F]), unlist(tmp[, ..expo, with=F]))
  
  ## run
  if(ncol(nphe) == 2 & nrow(nphe) == 3 & nphe[2,2] > 0){
    
    ## where are we at
    cat("run association between", expo, "and", outc, "\n")
    
    ## change adjustment
    adj <- ifelse(outc == "death.covid19.comprsk", gsub("pop", "non.european", adj), adj)
    
    ## run model: implement cox-prop test and possible step function for time intervals
    ff   <- tryCatch(
      {
        ## run multi-state cox model
        coxph(as.formula(paste0("Surv(", cdat,", as.factor(", outc, ")) ~ ", expo, " + ", adj)), data=tmp, ties = "breslow", id=f.eid)
      }, error=function(e){
        return(NA)
      })
    
    ## proceed only, if model converged
    if(length(ff) > 1){
      
      ## cox-prop test
      ff.p <- tryCatch(
        { ## test cox-proportional hazard assumption
          cox.zph(ff)$table[c(paste0(expo, "_1:2"), "GLOBAL"), 3]
        }, error=function(e){
          ## return NULL if untestable
          return(rep(NA, 2))
          
        })
      ## compute summary for storage
      ff   <- summary(ff)
      ## depend on outcome
      tmp  <- data.frame(outcome=outc, phecode=expo, beta=ff$coefficients[1,1], se=ff$coefficients[1,3], 
                         pval=ff$coefficients[1,6], nevent=sum(nphe[2,]), nall=sum(!is.na(tmp[, outc, with=F])), n.phecode=sum(nphe[,2]),
                         p.resid.phecode=ff.p[1], p.resid.overall=ff.p[2])
      ## return
      return(tmp)
    }
  }
}, mc.cores = 30) 
## combine into one
res.surv <- rbindlist(res.surv, fill = T)

## write results to file
fwrite(res.surv, paste0("output/Results.competing.", phe), sep = "\t", row.names = F, na = NA)
