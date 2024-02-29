#!/usr/bin/env Rscript

## script to perfrom MR for phecode on COVID-19 statistics
## Maik Pietzner 08/12/2022
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)
# print(R.Version())

setwd("<path to file>")

## --> packages needed <-- ##
require(MendelianRandomization)
require(mr.raps)
require(MRPRESSO)
require(data.table)
require(doMC)
# require(TwoSampleMR)

## --> import parameters <-- ##

## get regional coordinates
expo  <- args[1]

#-----------------------------------------#
##--       import relevant stats       --##
#-----------------------------------------#

## --> phecode statistics <-- ##

## import overall phecode results
res.phe   <- fread("<path to file>")
## subset to phecode of interest and omit MHC region
res.phe   <- res.phe[ id == expo & group != "6_11"]

## also include proxies
proxy.snp <- fread("<path to file>")
## keep only those relevant
proxy.snp <- proxy.snp[ MarkerName.1 %in% res.phe$MarkerName & R2.proxy >= .8]
## obtain those from the summary GWAS
res.tmp   <- fread(cmd=paste0("zcat <path to file>/", expo,".allchr.results.gz"))
## create MarkerName
res.tmp[, MarkerName := paste0("chr", CHROM, ":", GENPOS, "_", pmin(ALLELE0, ALLELE1), "_", pmax(ALLELE0, ALLELE1))]
## subset to those needed (add all stats of proxy SNPs, and choose later on based on largest LD)
res.tmp   <- merge(proxy.snp, res.tmp, by.x="MarkerName.2", by.y = "MarkerName")

## --> COVID-19 statistics <-- ##

## do in parallel
registerDoMC(3)

## import COVID-19 statistics (ALT - effect allele)
res.covid <- mclapply(c("A2", "B2", "C2"), function(x){
  
  ## import data
  tmp           <- fread(cmd=paste0("zcat gwas_stats/COVID19_HGI_", x, "_ALL_leave_23andme_and_UKBB_20220403_GRCh37.tsv.gz"))
  ## rename
  names(tmp)[1] <- "CHR"
  ## create MarkerName
  tmp[, MarkerName := paste0("chr", CHR, ":", POS, "_", pmin(REF, ALT), "_", pmax(REF, ALT))]
  ## subset to SNPs needed
  tmp           <- tmp[ MarkerName %in% res.tmp$MarkerName.2]
  ## add outcome
  tmp$outcome   <- x 
  ## return results
  return(tmp)
  
}, mc.cores = 3)
## combine into one large
res.covid <- do.call(rbind, res.covid)
## reshape
res.covid <- reshape(res.covid, idvar = c("MarkerName", "CHR", "POS", "REF", "ALT"), timevar = "outcome", direction = "wide")

#-----------------------------------------#
##--         combine and align         --##
#-----------------------------------------#

## combine both
res.comb <- merge(res.tmp, res.covid, by.x="MarkerName.2", by.y="MarkerName", suffixes=c(".phecode", ".covid"))
## decide on only one proxy per lead signal (MarkerName.1!)
res.comb <- res.comb[order(MarkerName.1, -R2.proxy)]
res.comb[, ind := 1:.N, by="MarkerName.1"]
res.comb <- res.comb[ ind == 1]

## keep only what is really needed
res.comb <- res.comb[, c("MarkerName.2", "ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P", "REF", "ALT", "R2.proxy",
                         grep("all_inv_var_meta_beta|all_inv_var_meta_sebeta|all_inv_var_meta_p|cases|controls", names(res.comb), value=T)), with=F]

## recode effects to be all phecode increasing
res.comb[, ea.phecode := ifelse(BETA > 0, ALLELE1, ALLELE0)]
res.comb[, nea.phecode := ifelse(BETA > 0, ALLELE0, ALLELE1)]
## recode beta accordingly
res.comb[, beta.phecode := abs(BETA)]
res.comb[, se.phecode := SE]

## convert to data frame to ease coding
res.comb <- as.data.frame(res.comb)

## align COVID-19 effects accordingly
for(j in c("A2", "B2", "C2")){
  res.comb[, paste0("beta.", j)] <- ifelse(res.comb$ALT == res.comb$ALLELE1, res.comb[, paste0("all_inv_var_meta_beta.", j)], -res.comb[, paste0("all_inv_var_meta_beta.", j)])
  res.comb[, paste0("se.", j)]   <- res.comb[, paste0("all_inv_var_meta_sebeta.", j)]
}

## keep only what is really needed
res.comb <- res.comb[, c("MarkerName.2", "ID", "ea.phecode", "nea.phecode", "beta.phecode", "se.phecode", "beta.A2", "se.A2",
                         "beta.B2", "se.B2", "beta.C2", "se.C2")]

## do some clean up
gc()

## some output along the way
print(res.comb)

## proceed only if at least five SNPs in the overlap
if(nrow(res.comb) >= 5){
  
  #-----------------------------------------#
  ##--           run MR methods          --##
  #-----------------------------------------#
  
  registerDoMC(3)
  
  ## loop through each phenotype
  res.mr <- mclapply(c("A2", "B2", "C2"), function(x){
    
    ## drop NAs
    tmp        <- na.omit(res.comb[, c("MarkerName.2", "beta.phecode", "se.phecode", paste0("beta.", x), paste0("se.", x))])
    
    ## --> MR-RAPS as first pass <-- ##
    
    ## run MR-PRAPS
    res.raps   <- mr.raps(b_exp = tmp$beta.phecode, b_out = tmp$se.phecode, 
                          se_exp = tmp[, paste0("beta.", x)], se_out = tmp[, paste0("se.", x)],
                          over.dispersion = T, loss.function = "huber")

    
    ## --> MR-PRESSO for pleiotropy outlier detection <-- ##
    
    ## MR-PRESSO
    res.presso <- mr_presso(BetaExposure = "beta.phecode", BetaOutcome = paste0("beta.", x),
                            SdExposure = "se.phecode", SdOutcome = paste0("se.", x), data=as.data.frame(tmp),
                            OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 10000, seed = 42)
    
    
    ## update SNP list if needed
    if(as.numeric(gsub("<", "", res.presso$`MR-PRESSO results`$`Global Test`$Pvalue)) < .05){
      
      ## test whether any outliers
      if(res.presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`[1] != "No significant outliers"){
        ## get the indices to be dropped
        ii        <- res.presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
        ## drop from tmp
        tmp       <- tmp[-ii, ]
      }
      
    }
    
    ## prepare input for MR models
    mr.input   <- mr_input(bx = tmp$beta.phecode, bxse = tmp$se.phecode,
                           by = tmp[, paste0("beta.", x)], byse = tmp[, paste0("se.", x)],
                           snps = tmp$MarkerName.2)
    
    ## run models
    res.mr     <- mr_allmethods(mr.input, method = "main")
    
    ## get estimate of heterogeneity
    res.het    <- mr_ivw(mr.input)@Heter.Stat
    
    #------------------------------------------#
    ##--  create common results data frame  --##
    #------------------------------------------#
    
    ## overall stats
    res.return                     <- data.frame(res.mr@Values)
    ## edit method
    res.return$Method              <- c("simple.median", "weighted.median", "ivw", "mr.egger", "mr.intercept")
    ## reshape
    res.return$expo                <- expo
    res.return                     <- reshape(res.return, idvar = "expo", timevar = "Method", direction = "wide")
    ## keep only what is of interest
    res.return                     <- res.return[, c("expo", grep("Estimate|Std.Err|P.value", names(res.return), value = T))]
    ## add heterogeneity
    res.return$het.i2              <- res.het[1]
    res.return$het.pval            <- res.het[2]
    
    ## add MR Raps
    res.return$Estimate.raps       <- res.raps$beta.hat 
    res.return$Std.Error.raps      <- res.raps$beta.se
    res.return$P.value.raps        <- res.raps$beta.p.value
    
    ## add MR Presso (all instruments)
    res.return$Estimate.presso     <- res.presso$`Main MR results`$`Causal Estimate`[1]
    res.return$Std.Error.presso    <- res.presso$`Main MR results`$Sd[1]
    res.return$P.value.presso      <- res.presso$`Main MR results`$`P-value`[1]
    ## add MR Presso (exclude instruments)
    res.return$Estimate.presso.wo  <- res.presso$`Main MR results`$`Causal Estimate`[2]
    res.return$Std.Error.presso.wo <- res.presso$`Main MR results`$Sd[2]
    res.return$P.value.presso.wo   <- res.presso$`Main MR results`$`P-value`[2]
    ## outlier test
    res.return$outlier.pval        <- res.presso$`MR-PRESSO results`$`Global Test`$Pvalue
    
    ## add phenotype
    res.return$pheno               <- x
    ## add number of SNPs
    res.return$n.snps              <- nrow(tmp)
    
    ## return results
    return(res.return)
    
  }, mc.cores=3)
  ## combine into one data frame
  res.mr <- do.call(rbind, res.mr)
  
  #-----------------------------------------#
  ##--             plot results          --##
  #-----------------------------------------#
  
  ## create confidence intervals to ease plotting
  res.comb$cil.phecode <- res.comb$beta.phecode - 1.96 * res.comb$se.phecode
  res.comb$ciu.phecode <- res.comb$beta.phecode + 1.96 * res.comb$se.phecode
  ## confidence intervals for COVID-19
  for(j in c("A2", "B2", "C2")){
    ## lower limit
    res.comb[, paste0("cil.", j)] <- res.comb[, paste0("beta.", j)] - 1.96 * res.comb[, paste0("se.", j)]
    ## upper limit
    res.comb[, paste0("ciu.", j)] <- res.comb[, paste0("beta.", j)] + 1.96 * res.comb[, paste0("se.", j)]
  }
  
  ## actual plot
  pdf(paste("graphics/MRresults", expo, "COVID19.pdf", sep="."), width = 6, height = 2)
  ## graphical parameters
  par(mar=c(1.5,1.5,.5,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, mfrow=c(1,3))
  
  ## loop over all three outcomes
  for(j in c("A2", "B2", "C2")){
    
    ## empty plot
    plot(range(res.comb[, c("cil.phecode", "ciu.phecode")]), range(res.comb[, paste0(c("cil.", "ciu."), j)], na.rm=T),
         type="n", ylab=paste0("Effect on ", j), xlab="Effect on Phecode")
    ## add zero crossings
    abline(v=0, lwd=.5); abline(h=0, lwd=.5)
    ## add confidence intervals
    arrows(res.comb$cil.phecode, res.comb[, paste0("beta.", j)], res.comb$ciu.phecode, res.comb[, paste0("beta.", j)], lwd=.4, length = 0,
           col="grey50")
    arrows(res.comb$beta.phecode, res.comb[, paste0("cil.", j)], res.comb$beta.phecode, res.comb[, paste0("ciu.", j)], lwd=.4, length = 0,
           col="grey50")
    ## add point estimates
    points(res.comb$beta.phecode, res.comb[, paste0("beta.", j)], pch=21, lwd=.3, bg="white", cex=.5)
    
    ## colour vector
    cl <- RColorBrewer::brewer.pal(6, "Set1")
    kk <- c("ivw", "raps", "presso", "simple.median", "weighted.median")
    ## add estimates
    for(k in 1:length(kk)){
      ## add the estimates
      abline(a=0, b=res.mr[which(res.mr$pheno == j), paste0("Estimate.", kk[k])], lwd=.5, col=cl[k])
    }
    ## add MR Egger
    abline(a=res.mr$Estimate.mr.intercept[which(res.mr$pheno == j)], 
           b=res.mr$Estimate.mr.egger[which(res.mr$pheno == j)], lwd=.5, col=cl[6])
    
    ## add legend
    legend("bottomright", bty="n", cex=.3, pch=NA, col=cl, 
           legend = c("ivw", "raps", "presso", "simple.median", "weighted.median", "mr.egger"),
           lty=1)
  }
  
  dev.off()
  
  #------------------------------------------#
  ##--            export results          --##
  #------------------------------------------#
  
  write.table(res.mr, paste("output/MRresults", expo, "COVID19", sep="."), sep="\t", row.names = F)
  
}
