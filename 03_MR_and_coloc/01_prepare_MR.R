######################################################
#### MR phecodes --> COVID-19                     ####
#### Maik Pietzner                                ####
######################################################

rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)
load(".RData")

## packages needed
require(data.table)
# require(readxl)


#######################################################
####      import genetic correlation estimates     ####
#######################################################

## import phecode labels
phecode.label <- fread("<path to file>")

## import genetic correlation done by Summaria
res.cor       <- fread("<path to file>/covid_release7_vs_allphecodes.LF.20221214")
## use same poplulation as for the MR
res.cor       <- res.cor[ p1 %in% c("A2_ALL_exUKBB", "B2_ALL_exUKBB", "C2_ALL_exUKBB") & p2 %in% phecode.label$id & !is.na(rg)]
## select only those at least three times in
ii            <- table(res.cor$p2)
ii            <- names(ii[ ii == 3])
## subset again
res.cor       <- res.cor[ p2 %in% ii]
phecode.label <- phecode.label[ id %in% ii]

## add q-value
res.cor[, qval := p.adjust(p, method = "BH")]
## add label
res.cor       <- merge(res.cor, phecode.label[, c("id", "phecode", "phenotype", "category", "cl", "y.pos", "cases", "controls", "case.ratio.pcare")], by.x="p2", by.y="id")

#######################################################
####   create input file for phecode --> COVID-19  ####
#######################################################

## at least 5 sentinel variants
write.table(phecode.label[ number.sentinel >= 5]$id, "mr.input.phecode.covid19.txt", sep="\t", row.names = F, col.names = F, quote = F)

######################################
####   import reverse MR results  ####
######################################

## import output
ii         <- dir("../output/")
## subset to those related to phecodes as outcome
ii         <- grep("COVID19", ii, value=T)

## import results
res.rev.mr <- lapply(ii, function(x){
  
  ## import 
  tmp       <- fread(paste0("../output/", x))
  return(tmp)
  
})
## combine into one file
res.rev.mr <- do.call(rbind, res.rev.mr)
## add label
res.rev.mr <- merge(phecode.label[, c("id", "phecode", "phenotype", "category", "cl", "y.pos", "cases", "controls")], res.rev.mr, 
                    by.x="id", by.y="expo")

#-----------------------------------#
##--    prepare for manuscript   --## 
#-----------------------------------#

## make mr-presso p-value criterion to distinguish
res.rev.mr[, P.value.presso.combined := ifelse(!is.na(P.value.presso.wo), P.value.presso.wo, P.value.presso)]
## add multiple testing correction
res.rev.mr[, Q.value.presso.combined := p.adjust(P.value.presso.combined, method = "BH")]
## how many findings
nrow(res.rev.mr[ Q.value.presso.combined < .05])
## none...; apart from possibly AD

## write results to file
write.table(res.rev.mr, "<path to file>", sep="\t", row.names=F)
