############################################################
#### Cox-proportional hazard models for COVID-19        ####
#### Maik Pietzner                           31/08/2023 ####
############################################################

rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(arrow)
library(tidyr)
library(dplyr)
require(doMC)
require(survival)
require(basicPlotteR)
require(colorspace)
require(igraph)
require(tidyverse)
require(poweRlaw)
require(metafor)

##################################################
####         import basic covariates          ####
##################################################

## An example list of columns:
cl.select       <- c("f.eid", "f.21022.0.0", "f.31.0.0", "f.53.0.0", "f.54.0.0", "f.21001.0.0", paste0("f.22009.0.", 1:10))
## names to assign
cl.names        <- c("f.eid", "age", "sex", "baseline_date", "centre", "bmi", paste0("pc", 1:10))

## import data from the main release
ukb.dat         <- read_parquet("<path to file>",
                                col_select = cl.select)
## change names
names(ukb.dat)  <- cl.names

## import withdrawals
withdrawls      <- fread("<path to file>")
ukb.dat         <- subset(ukb.dat, !(f.eid %in% withdrawls$V1))

## import ancestry information (loop over all)
ukb.anc         <- lapply(c("EUR", "AFR", "CSA", "EAS", "MID", "AMR"), function(x){
  ## import assignment
  tmp <- fread(paste0("<path to file>", x, "_panukbb.ids"))
})
## combine all into one large file
ukb.anc         <- do.call(rbind, ukb.anc)
## streamline population assignment
ukb.anc[, pop := gsub("NA_", "", pop)]

##################################################
####                death data                ####
##################################################

## import death dates for censoring
ukb.dod  <- fread("<path to file>")
## add date of death (keep only first instance)
ukb.comb <- merge(ukb.comb, ukb.dod[ ins_index == 0 , c("eid", "date_of_death")], all.x=T, by.x="f.eid", by.y="eid")

## how many died before the pandemic
nrow(ukb.dod[ as.Date(date_of_death, format="%d/%m/%Y") < as.Date("01/01/2020", format = "%d/%m/%Y")])
## n = 28,606

##################################################
####          COVID-19 testing data           ####
##################################################

## --> with support from Tomoko Nakanishi <-- ##

## import testing results
covid1   <- fread("<path to file>/covid19_result_england_RAW_230829.csv")
covid2   <- fread("<path to file>/covid19_result_scotland_RAW_230829.csv")
covid3   <- fread("<path to file>/covid19_result_wales_RAW_230829.csv")

## simplify to include only cases and earliest date
covid1   <- covid1 %>% filter(result == 1) %>% group_by(eid) %>%
  mutate(hospital = max(hosaq),
         covid_date = min(as.Date(specdate, format="%d/%m/%Y"))) %>%
  distinct(eid, .keep_all=TRUE)

## further simplify
covid1   <- covid1 %>% mutate(hospital = ifelse(hospital==1, 1, 0))
covid1   <- covid1 %>% select(eid, result, hospital, covid_date)
## n = 109,404
summary(covid1)

## same for Scottish data
covid2   <- covid2 %>% filter(result == 1) %>% group_by(eid) %>%
  mutate(hospital = ifelse(any(factype == 3), 1, 0),
         covid_date = min(as.Date(specdate, format="%d/%m/%Y"))) %>%
  distinct(eid, .keep_all=TRUE)
covid2   <- covid2 %>% select(eid, result, hospital, covid_date)
## n = 1,679
summary(covid2)

## and Welsh data
covid3   <- covid3 %>% filter(result == 1) %>% group_by(eid) %>%
  mutate(hospital = ifelse(any(pattype == 7 & !(perstype %in% c(1:38,76:77,82:94,100:101))), 1, 0),
         covid_date = min(as.Date(specdate, format="%d/%m/%Y"))) %>%
  distinct(eid, .keep_all=TRUE)
covid3   <- covid3 %>% select(eid, result, hospital, covid_date)
## n = 3,171
summary(covid3)

## combine
covid    <- bind_rows(covid1, covid2, covid3)
covid    <- covid %>% arrange(covid_date)
## n = 114,254

## clean to avoid duplications
covid    <- covid %>% group_by(eid) %>%
  mutate(hospital = max(hospital),
         covid_date = min(covid_date)) %>%
  distinct(eid, .keep_all=TRUE)
## n = 111,601

## convert 
covid    <- as.data.table(covid)

## add to the data
ukb.comb <- merge(ukb.comb, covid, all.x = T, by.x="f.eid", by.y="eid")
ukb.comb <- as.data.table(ukb.comb)

## creaete addiotinal variables
ukb.comb[, age_at_diagnosis := age + round(as.numeric(as.Date(covid_date) - as.Date(baseline_date))/365)]
summary(ukb.comb$age_at_diagnosis)
## still living at the end of the study
ukb.comb[, death := ifelse(is.na(date_of_death), 0, 1)]
## whether reported COVID-19
ukb.comb[, sarscov2.status := ifelse(is.na(covid_date), 0, 1)]
table(ukb.comb$sarscov2.status)
#      0      1 
# 362052 103397 

##################################################
####    Hospitalisation/Death/LongCOVID       ####
##################################################

## import HES data
hesin_diag    <- fread("<path to file>/hesin_diag.txt", sep="\t")
hesin         <- fread("<path to file>/hesin.txt", sep="\t")

## look at participants tested negative
neg           <- unique(ukb.comb[sarscov2.status == 0]$f.eid)
tmp           <- hesin[hesin$eid %in% neg,]
## merge visit with diagnosis
hesin_new_neg <- merge(tmp, hesin_diag, by=c("eid", "ins_index"))
## filter for hospitalisation with COVID-19
hesin_new_neg <- hesin_new_neg %>% filter(diag_icd10 %in% c("U071","U072"))
## order by date, to get earliest occurrence
hesin_new_neg <- hesin_new_neg[order(as.Date(epistart, format="%d/%m/%Y")),]
## drop repeated hospitalisations
hesin_new_neg <- hesin_new_neg[!duplicated(eid),]

## recode hospital cases in the data
for(i in seq(1, length(hesin_new_neg$eid))){
  ukb.comb$sarscov2.status[ukb.comb$f.eid == hesin_new_neg$eid[i]] <- 1
  ukb.comb$hospital[ukb.comb$f.eid == hesin_new_neg$eid[i]]        <- 1
  ukb.comb$covid_date[ukb.comb$f.eid == hesin_new_neg$eid[i]]      <- as.Date(hesin_new_neg$epistart[i], format="%d/%m/%Y")
}


## look at people with positive Sars-CoV2 status
hesin_new     <- merge(hesin, hesin_diag, by=c("eid", "ins_index"))
## only events after the pandemic hit
hesin_new     <- hesin_new %>% filter(as.Date(epistart, format="%d/%m/%Y") >= as.Date("2020-01-01"))
## filter to people with positive tests
hesin_new     <- hesin_new %>% filter(eid %in% ukb.comb[!is.na(covid_date)]$f.eid)
## reduce to COVID-19 codes
hesin_covid   <- hesin_new %>% filter(diag_icd10 %in% c("U071", "U072"))

## look at procedures after the pandemic
hesin_oper    <- fread("<path to file>hesin_oper.txt", sep="\t")
hesin_oper    <- hesin_oper %>% filter(as.Date(opdate, format="%d/%m/%Y") >= as.Date("2020-01-01"))
## combine with HES data
hesin_new1    <- merge(hesin, hesin_oper, by=c("eid", "ins_index"))
hesin_new1    <- hesin_new1 %>% filter(eid %in% ukb.comb[!is.na(covid_date)]$f.eid)

## critical care admitted to ICU within -14 and 30 from the COVID-19 diagnosis date
c             <- fread("<path to file>/hesin_critical.txt", sep="\t")
c             <- c[c$eid %in% ukb.comb$f.eid,]
## add data
c2            <- merge(ukb.comb, c, by.x = "f.eid", by.y = "eid", all.x = T)
## filter
c2            <- c2 %>% filter(as.Date(ccstartdate, "%d/%m/%Y") <= 30 + covid_date)
c2            <- c2 %>% filter(as.Date(ccstartdate, "%d/%m/%Y") >= covid_date - 14)
c2            <- c2 %>% filter(aressupdays > 0)
## get the relevant IDs for those requiring severe support
c2id          <- unique(c2$f.eid)

## add column to the data
ukb.comb[, severe_resp := ifelse(f.eid %in% c2id & sarscov2.status == 1, 1, 0)]
table(ukb.comb$severe_resp)

## use J80, J9600, J9609, Z991 in HES for the severe respiratory failure within -14 and 30 from the COVID-19 diagnosis date
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- ukb.comb[ukb.comb$f.eid == unique(hesin_new$eid)[i],]
  if(dim(TMP)[1]>0){
    tmp1 <- tmp %>% filter(as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid_date) + 30)
    tmp1 <- tmp1 %>% filter(as.Date(epistart, "%d/%m/%Y") >= as.Date(TMP$covid_date) - 14 | as.Date(admidate, "%d/%m/%Y") >= as.Date(TMP$covid_date) - 14)
    tmp1 <- tmp1 %>% filter(as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid_date))
    ukb.comb$severe_resp[ukb.comb$f.eid == unique(hesin_new$eid)[i] & "J80" %in% tmp1$diag_icd10]   <- 1
    ukb.comb$severe_resp[ukb.comb$f.eid == unique(hesin_new$eid)[i] & "J9600" %in% tmp1$diag_icd10] <- 1
    ukb.comb$severe_resp[ukb.comb$f.eid == unique(hesin_new$eid)[i] & "J9609" %in% tmp1$diag_icd10] <- 1
    ukb.comb$severe_resp[ukb.comb$f.eid == unique(hesin_new$eid)[i] & "Z991" %in% tmp1$diag_icd10]  <- 1
  }
}

#use E851, E852 in OPER code for the severe respiratory failure within -14 and 30 from the COVID-19 diagnosis date
for(i in seq(1,length(unique(hesin_new1$eid)))){
  tmp <- hesin_new1[hesin_new1$eid == unique(hesin_new1$eid)[i],]
  TMP <- ukb.comb[ukb.comb$f.eid == unique(hesin_new1$eid)[i],]
  if(dim(TMP)[1]>0){
    tmp1 <- tmp %>% filter(as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid_date) + 30)
    tmp1 <- tmp1 %>% filter(as.Date(epistart, "%d/%m/%Y") >= as.Date(TMP$covid_date) - 14 | as.Date(admidate, "%d/%m/%Y") >= as.Date(TMP$covid_date) - 14)
    tmp1 <- tmp1 %>%  filter(as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid_date))
    ukb.comb$severe_resp[ukb.comb$f.eid == unique(hesin_new1$eid)[i] & "E851" %in% tmp1$oper4] <- 1
    ukb.comb$severe_resp[ukb.comb$f.eid == unique(hesin_new1$eid)[i] & "E852" %in% tmp1$oper4] <- 1
  }
}

## --> define time data <-- ##

## define right censoring
ukb.comb[, covid.date.analysis := ifelse(!is.na(covid_date), as.character(covid_date), "2022-12-31")]
ukb.comb[, covid.date.analysis := as.Date(covid.date.analysis)]
summary(ukb.comb$covid.date.analysis)

## do some cleaning
rm(list=c(grep("hesin*", ls(), value=T), grep("^covid*", ls(), value = T)))
gc(reset=T)

## --> COVID-19 related death <-- ##

## import death certificates
ukb.dead <- fread("<path to file>/death_cause.txt")
## subset to those with evidence of COVID-19 on the certificate
ukb.dead <- ukb.dead[ cause_icd10 %in% c("U071", "U072")]
length(unique(ukb.dead$eid))

## add to the data
ukb.comb <- merge(ukb.comb, ukb.dead, by.x="f.eid", by.y="eid", all.x=T)
## create new variable
ukb.comb[, death.covid19 := ifelse(!is.na(cause_icd10), 1, 0)]
## drop some
ukb.comb$cause_icd10 <- ukb.comb$ins_index <- ukb.comb$arr_index <- ukb.comb$level <- NULL
## define death censoring
summary(as.Date(ukb.dod$date_of_death, "%d/%m/%Y"))
ukb.comb[, date.death.covid19 := ifelse(!is.na(date_of_death), date_of_death, "31/12/2022")]
ukb.comb[, date.death.covid19 := as.Date(date.death.covid19, "%d/%m/%Y")]
summary(ukb.comb$date.death.covid19)

## --> Long COVID (primary care data) <-- ##

## codes to be used
long.covid        <- c("Y2b89", "Y2b8a", "Y2b87", "Y2b88", "1325161000000102", "1325031000000108", "1325041000000104", "1325181000000106", "1325021000000106",
                       "1325141000000103", "1325081000000107", "1325061000000103", "1325071000000105", "1325051000000101") 

## import TPP data
ukb.tpp           <- fread("<path to file>/covid19_tpp_gp_clinical.txt")
## reduce to the set of relevant codes: CTV3 codes: "Y2b89" "Y2b8a" "Y2b87" "Y2b88"
ukb.tpp           <- ukb.tpp[ code %in% long.covid]
## reduce to earliest diagnosis
ukb.tpp           <- ukb.tpp[ order(eid, event_dt)]
ukb.tpp[, ind := 1:.N, by="eid"]
ukb.tpp           <- ukb.tpp[ ind == 1 ]
summary(ukb.tpp)
gc(reset=T)

## import EMIS data
ukb.emis          <- fread("<path to file>/covid19_emis_gp_clinical.txt")
## reduce to relevant codes: SNOMED CT: "1325161000000102" "1325031000000108" "1325041000000104" "1325181000000106" "1325021000000106" "1325141000000103" "1325081000000107" "1325061000000103" "1325071000000105" "1325051000000101"
ukb.emis          <- ukb.emis[ code %in% long.covid]
## reduce to earliest diagnosis
ukb.emis          <- ukb.emis[ order(eid, event_dt)]
ukb.emis[, ind := 1:.N, by="eid"]
ukb.emis          <- ukb.emis[ ind == 1 ]

## combine
long.covid        <- unique(rbind(ukb.tpp[, c("eid", "event_dt", "ind")], ukb.emis[, c("eid", "event_dt", "ind")]))
length(unique(long.covid$eid))

## rename
names(long.covid) <- c("f.eid", "long.covid.date", "long.covid.status")

## add to the general data
ukb.comb          <- merge(ukb.comb, long.covid, all.x=T)
## recode
ukb.comb[, long.covid.status := ifelse(is.na(long.covid.status), 0, 1)]
## needs probably an earlier ending...
ukb.comb[, long.covid.date := ifelse(is.na(long.covid.date), "30/09/2021", long.covid.date)]
ukb.comb[, long.covid.date := as.Date(long.covid.date, format = "%d/%m/%Y")]
summary(ukb.comb$long.covid.date)
## delete one person with implausible date from analysis
ukb.comb[ long.covid.date < as.Date("01/01/2020", format = "%d/%m/%Y")]$long.covid.status <- NA

## --> additional definition of hospitalisation <-- ##

## import HES data
hesin_diag    <- fread("<path to file>/hes_data_20230406/hesin_diag.txt", sep="\t")
hesin         <- fread("<path to file>/hes_data_20230406/hesin.txt", sep="\t")
## observation period
summary(as.Date(hesin$epistart, format = "%d/%m/%Y"))

## subset to relevant codes
hesin_diag    <- hesin_diag[ diag_icd10 %in% c("U071", "U072") ]
## add relevant times
hesin_diag    <- merge(hesin, hesin_diag, by=c("eid", "ins_index"))
## drop dates before the pandemic hit
hesin_diag    <- hesin_diag[ as.Date(epistart, format="%d/%m/%Y") >= as.Date("2020-01-01") ]
## subset to first occurence
hesin_diag    <- hesin_diag[order(as.Date(epistart, format="%d/%m/%Y"))]
hesin_diag    <- hesin_diag[ !duplicated(eid) ]

## add to the data
ukb.comb      <- merge(ukb.comb, hesin_diag[, c("eid", "epistart")], by.x="f.eid", by.y="eid", all.x=T)
## create new variable
ukb.comb[, hosp.covid19 := ifelse(!is.na(epistart), 1, 0)]
## define censoring
ukb.comb[, date.hosp.covid19 := ifelse(!is.na(epistart), epistart, "31/12/2022")]
ukb.comb[, date.hosp.covid19 := as.Date(date.hosp.covid19, "%d/%m/%Y")]
table(ukb.comb$hosp.covid19)

##################################################
####          general data tweaking           ####
##################################################

## delete what is no longer needed
rm(list=setdiff(ls(), "ukb.comb"))
gc(reset=T)

## drop participants never exposed to the virus
ukb.comb <- ukb.comb[ is.na(date_of_death) | date_of_death > as.Date("01/03/2020")]

## create artificial baseline date to define prevalent health records up to here
ukb.comb[, baseline.covid := as.Date("01/01/2020", format = "%d/%m/%Y")]

## age as the pandemic hit
ukb.comb[, age.covid := age + round(as.numeric(baseline.covid - as.Date(baseline_date))/365) ]
summary(ukb.comb$age.covid)

##################################################
####            import phecodes               ####
##################################################

## import phecode data
require(data.table)
ukbb.phe     <- fread("<path to file>/UKBB.phecodes.collated.sources.20210811.txt", sep="\t", header = T)

## convert to wide format
ukbb.phe     <- dcast(ukbb.phe, f.eid ~ phecode, value.var = c("date", "resource"))

## import label
lab.phe      <- read.table("<path to file>Comparison.phecodes.cases.UKBB.UCL.Charite.20210811.txt", sep="\t", header=T)
## add identifier to match with data set
lab.phe$id   <- paste0("date_", lab.phe$phecode) 
## drop all with less than 100 cases
lab.phe      <- subset(lab.phe, n.charite >= 100) 
## N = 1,448

## add sex-specific coding
lab.tmp      <- read.csv("<path to file>/phecode_definitions1.2.csv")
## create variable to map to UKBB data set
lab.tmp$id   <- paste0("date_", lab.tmp$phecode) 
## add sex-specific information to the data
lab.phe      <- merge(lab.phe, lab.tmp[, c("id", "sex")], by="id")
## replace missing ones
lab.phe$sex[lab.phe$sex == ""] <- "Both"

## restrict to what is really needed
ukbb.phe     <- ukbb.phe[, c("f.eid", lab.phe$id), with=F]

## add to covariates (also to look for further controls)
ukb.comb     <- merge(ukb.comb, ukbb.phe, all.x=T, by="f.eid")

## convert dates to binary variables based on comparison with baseline date
ukb.comb               <- as.data.frame(ukb.comb)
ukb.comb[, lab.phe$id] <- apply(ukb.comb[, lab.phe$id], 2, function(x){
  x <- ifelse(!is.na(x) & as.Date(x) < as.Date("2020-03-01"), 1, 0)
  return(x)
})

## clean up
rm(ukbb.phe); rm(lab.tmp)

##################################################
####    association testing w/ primary care   ####
##################################################

## convert to data table
ukb.comb <- as.data.table(ukb.comb)

## compute follow-up time for each event (01/01/2020 as start date); severe resp
ukb.comb[, fol.covid := as.numeric((covid.date.analysis - baseline.covid))/365.25]
summary(ukb.comb$fol.covid)

## hospitalisation
ukb.comb[, fol.hosp := as.numeric((date.hosp.covid19 - baseline.covid))/365.25]
summary(ukb.comb$fol.hosp)

## death
ukb.comb[, fol.death := as.numeric((date.death.covid19 - baseline.covid))/365.25]
summary(ukb.comb$fol.death)
## drop people died before the pandemic
ukb.comb <- ukb.comb[ fol.death > 0]

## date long covid
ukb.comb[, fol.long.covid := as.numeric((long.covid.date - baseline.covid))/365.25]

## simplify sex for analysis
ukb.comb[, sex.bin := ifelse(sex == "Female", 0, 1)]
## define EUR as a reference population
ukb.comb[, pop := factor(pop, levels = c("EUR", "AFR", "AMR", "CSA", "EAS", "MID"))]

#----------------------------------#
##--        across diseases     --##
#----------------------------------#

## write to file for association testing
fwrite(lab.phe, "Labels.phecodes.20230904.txt", sep="\t", na = NA, row.names = F)
fwrite(ukb.comb, "UKB.phecodes.COVID19.20230904.txt", sep="\t", na = NA, row.names = F)

#----------------------------------#
##--          sum numbers       --##
#----------------------------------#

## follow-up time summary
ukb.comb <- as.data.table(ukb.comb)
summary(ukb.comb[ hosp.covid19 == 1]$date.hosp.covid19)

##################################################
####             import results               ####
##################################################

#------------------------#
##-- w/  primary care --##
#------------------------#

## import
res.w.primary  <- fread("../output/Results.UKB.phecodes.COVID19.20230904.txt")
## add label
res.w.primary  <- merge(lab.phe, res.w.primary, by.x = "id", by.y = "phecode")
res.w.primary  <- as.data.table(res.w.primary)


## drop two large data sets to save memory
rm(ukb.comb); rm(ukb.comb2); gc(reset=T)

## subset labels
lab.cox        <- lab.phe[ id %in% res.w.primary$id]
## sort
lab.cox        <- lab.cox[order(phecode)]
## add new variable to sort in the plot latter on
lab.cox[, srt := 1:nrow(lab.cox)]
## add to the results
res.w.primary  <- merge(res.w.primary, lab.cox[, c("id", "srt")], by="id")
## convert to data table
res.w.primary  <- as.data.table(res.w.primary)

## number of significant findings
nrow(subset(res.w.primary, pval < .05/nrow(res.w.primary))) ## N = 1148

## number of diseases across all outcomes
length(unique(subset(res.w.primary, pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3 )$phecode)) ## N = 679

## only severe COVID-19
length(unique(subset(res.w.primary, pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3 & outcome != "long.covid.status")$phecode)) ## N = 679

#------------------------#
##--     figures      --##
#------------------------#

## create plotting vector
cl.cat        <- do.call(data.frame, aggregate(srt ~ category, lab.cox, function(x) c(min(x), mean(x), max(x))))
names(cl.cat) <- c("category", "start", "mid", "end")
## add colour coding
cl.cat        <- merge(cl.cat, unique(lab.cox[, c("category", "cl")]))
## order
cl.cat        <- cl.cat[ order(cl.cat$start), ]

## create label for COVID-19 outcomes
lab.covid     <- data.frame(covid=c("hosp.covid19", "severe_resp", "death.covid19", "long.covid.status"), 
                            date=c("fol.hosp", "fol.covid", "fol.death", "fol.long.covid"),
                            label=c("Hospitalisation", "Respiratory failure", "Death", "LongCOVID"),
                            srt=1:4)

## --> w/ primary care data <-- ##

## open device
pdf("../graphics/Summary.phecodes.COVID19.w.primary.20230907.pdf", width = 6.3, height = 3.94)
## plotting coordinates
par(mar=c(.2,1.5,.5,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
## layout
layout(matrix(1:4,4,1), heights = c(.2,.2,.2,.3))

## run through each phenotype
for( j in 1:nrow(lab.covid)){
  
  ## adopt plotting range
  if(j == 4) par(mar=c(3.5,1.5,.5,.5))
  
  ## empty plot
  plot(c(.5, nrow(lab.cox)+.5), c(0, max(-log10(res.combined[ outcome == lab.covid$covid[j]]$pval.w.primary))*1.1),
       type = "n", xlab="", xaxt="n", ylab=expression(log[10]("p-value")), yaxt="n")
  ## add axis
  axis(2, lwd=.5)
  ## get plotting coordinates to ease downstream coding
  pm  <- par("usr")
  ## add coloured background
  rect(cl.cat$start-.5, pm[3], cl.cat$end+.5, pm[4], border = NA, col = lighten(cl.cat$cl, .6))
  
  ## p-value threshold
  abline(h=-log10(.05/nrow(res.combined)), lwd=.5, lty=2)
  
  ## add associations
  tmp <- res.combined[ outcome == lab.covid$covid[j] ]
  points(tmp$srt, -log10(tmp$pval.w.primary), 
         cex=ifelse(tmp$pval.w.primary < .05/nrow(res.combined), .5, .3),
         pch=ifelse(tmp$beta.w.primary > 0, 24, 25),
         col="white",
         lwd=.2, bg="grey50")
  
  ## add disease label
  text(pm[1], pm[4]-(pm[4]-pm[3])*.05, labels = lab.covid$label[j], pos=4, cex=.5)
  
  ## top five annotations
  tmp <- tmp[ order( pval.w.primary)]
  tmp <- tmp[1:10]
  ## add to the plot
  addTextLabels(tmp$srt, -log10(tmp$pval.w.primary), tmp$phenotype, cex.label = .5,
                col.label="grey20", lwd = .2, keepLabelsInside = T, cex.pt = .7)
  
  ## add legend
  if(j == 4){
    ## axis label
    mtext("Diseases ordered by ICD-10 chapter", 1, cex=.6, line=.5)
    ## legend for phecode categories
    ## add legend
    legend(pm[1]-(pm[2]-pm[1])*.02, pm[3]-(pm[4]-pm[3])*.25, lty=0, pch=22, pt.lwd=.1, bty="n", xpd=NA,
           pt.bg=cl.cat$cl, legend = c(stringr::str_to_sentence(cl.cat$category[-nrow(cl.cat)]), "Other"),
           cex=.5, ncol=9, pt.cex=1)
  } 
}
## close device
dev.off()

#------------------------#
##--      numbers     --##
#------------------------#

## number of significant associations
nrow(res.w.primary[ pval < .05/nrow(res.w.primary)]) ## n = 1148
## how many passing proportional hazard assumption
nrow(res.w.primary[ pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3]) ## 1128

## add hazard ratio column
res.w.primary[, hr := exp(beta)]

## prepare for table
res.w.primary[, hr_column := paste0(sprintf("%.2f", hr), " (", sprintf("%.2f", exp(beta - 1.96*se)),";", sprintf("%.2f", exp(beta + 1.96*se)),")")]

## create a file in wide format
res.wide  <- reshape(res.w.primary, idvar = names(res.w.primary)[1:11], timevar = "outcome", direction = "wide")
## n = 1396

## compute meta-analysis across severe outcomes
res.wide  <- lapply(1:nrow(res.wide), function(x){
  
  ## compute MA across COVID-19 outcomes
  ma.cov <- rma(yi=c(res.wide$beta.death.covid19[x], res.wide$beta.hosp.covid19[x], res.wide$beta.severe_resp[x]),
                sei=c(res.wide$se.death.covid19[x], res.wide$se.hosp.covid19[x], res.wide$se.severe_resp[x]),
                method="FE")
  
  ## return results
  return(data.table(res.wide[x, ], heterogeneity=ma.cov$I2, pval.hetero=ma.cov$QEp))
  
})
## combine again
res.wide  <- do.call(rbind, res.wide)
res.wide  <- as.data.table(res.wide)

## --> look at diseases with heterogeneity <-- ##

## how many diseases with little evidence of hetero
nrow(res.wide[ ((pval.hosp.covid19 < .05/nrow(res.w.primary) & p.resid.phecode.hosp.covid19 > 1e-3) |
                  (pval.death.covid19 < .05/nrow(res.w.primary) & p.resid.phecode.death.covid19 > 1e-3) |
                  (pval.severe_resp < .05/nrow(res.w.primary) & p.resid.phecode.severe_resp > 1e-3)) &
                 pval.hetero > 1e-3])
nrow(res.wide[ ((pval.hosp.covid19 < .05/nrow(res.w.primary) & p.resid.phecode.hosp.covid19 > 1e-3) |
                  (pval.death.covid19 < .05/nrow(res.w.primary) & p.resid.phecode.death.covid19 > 1e-3) |
                  (pval.severe_resp < .05/nrow(res.w.primary) & p.resid.phecode.severe_resp > 1e-3))])
## n = 641 out of 672

## --> write to file for Supplemental table 3 <-- ##
write.table(res.wide, "Results.phecodes.COVID19.outcome.20230914.txt", sep="\t", row.names = F)

## enrichment of certain categories
res.enric <- as.data.frame(expand.grid(covid=lab.covid$covid, category=cl.cat$category))
res.enric <- lapply(1:nrow(res.enric), function(x){
  
  ## get outcome
  outc <- res.enric$covid[x]
  ## category
  dcat <- res.enric$category[x]
  
  ## sig diseases in category
  d1  <- nrow(res.w.primary[ outcome == outc & category == dcat & pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3])
  ## non-sig diseases in category
  d2  <- nrow(res.w.primary[ outcome == outc & category == dcat & pval > .05/nrow(res.w.primary) & p.resid.phecode > 1e-3])
  ## sig diseases not in category
  d3  <- nrow(res.w.primary[ outcome == outc & category != dcat & pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3])
  ## non-sig diseases not in category
  d4  <- nrow(res.w.primary[ outcome == outc & category != dcat & pval > .05/nrow(res.w.primary) & p.resid.phecode > 1e-3])
  
  ## test for enrichment
  enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
  
  ## return information needed
  return(data.frame(res.enric[x,], or=enr$estimate, pval=enr$p.value))
  
})
## combine
res.enric <- do.call(rbind, res.enric)

##################################################
####    Effect modification by ancestry/sex   ####
##################################################

## import results for possible effect modification
res.inter <- lapply(c("sex", "non.european"), function(x){
  
  ## import results
  tmp <- fread(paste0("../output/Results.", x,".UKB.phecodes.COVID19.20230904.txt"))
  ## add interaction type
  tmp[, term := x]
  ## return
  return(tmp)
  
})
## combine
res.inter <- do.call(rbind, res.inter)
## add label
res.inter <- merge(lab.phe, res.inter, by.x = "id", by.y = "phecode")
res.inter <- as.data.table(res.inter)

## significant findings
res.inter[ pval.inter < .05/nrow(res.inter) ]
## n = 9

## sanity check 
cox.tmp <- coxph(Surv(fol.hosp, hosp.covid19) ~  date_274.1*sex + age.covid + sex + non.european, ukb.comb, ties = "breslow")
summary(cox.tmp)
cox.tmp <- coxph(Surv(fol.hosp, hosp.covid19) ~  date_274.1 + age.covid + non.european, subset(ukb.comb, sex == "Male"), ties = "breslow")
summary(cox.tmp)
cox.tmp <- coxph(Surv(fol.hosp, hosp.covid19) ~  date_274.1 + age.covid + non.european, subset(ukb.comb, sex != "Male"), ties = "breslow")
summary(cox.tmp)

## create HR columns
res.inter[, hr_column.1 := paste0(sprintf("%.2f", exp(beta.1)), " (", sprintf("%.2f", exp(beta.1 - 1.96*se.1)),";", sprintf("%.2f", exp(beta.1 + 1.96*se.1)),")")]
res.inter[, hr_column.2 := paste0(sprintf("%.2f", exp(beta.2)), " (", sprintf("%.2f", exp(beta.2 - 1.96*se.2)),";", sprintf("%.2f", exp(beta.2 + 1.96*se.2)),")")]

## write to file
write.table(res.inter[ pval.inter < .05/nrow(res.inter)], "Sig.Results.phecode.COVID19.sex.euro.interaction.20230929.txt", sep="\t", row.names = F)

##################################################
####       disease - disease interaction      ####
##################################################

## write to file
write.table(expand.grid("UKB.phecodes.COVID19.20230904.txt", lab.cox$id), 
            "inpute.phecode.interaction.txt", sep="\t", row.names = F, col.names = F, quote = F)

## get output
ii <- dir("../output/")
ii <- grep("Results.d", ii, value=T)

## import results for possible effect modification
res.inter.phecode <- mclapply(ii, function(x){
  
  ## import results
  tmp <- fread(paste0("../output/", x))
  ## add interaction type
  tmp[, term := gsub("Results\\.|\\.UKB.phecodes.COVID19.20230904.txt", "", x)]
  ## return
  return(tmp)
  
}, mc.cores=10)
## combine
res.inter.phecode <- do.call(rbind, res.inter.phecode)
## add label (for primary exposure)
res.inter.phecode <- merge(res.inter.phecode, lab.phe[, c("id", "phenotype", "category")], by.x = "phecode", by.y = "id")
## add label (for interaction)
res.inter.phecode <- merge(res.inter.phecode, lab.phe[, c("id", "phenotype", "category")], by.x = "term", by.y = "id", suffixes = c(".primary", ".inter"))
## convert to data table
res.inter.phecode <- as.data.table(res.inter.phecode)

## sanity check since some findings look odd
cox.tmp <- coxph(Surv(fol.hosp, hosp.covid19) ~  date_709.7*date_1010 + age.covid + sex + non.european, ukb.comb, ties = "breslow")
summary(cox.tmp)
cox.tmp <- coxph(Surv(fol.hosp, hosp.covid19) ~  date_709.7 + age.covid + sex + non.european, subset(ukb.comb, date_1010 == 0), ties = "breslow")
summary(cox.tmp)
cox.tmp <- coxph(Surv(fol.hosp, hosp.covid19) ~  date_709.7 + age.covid + sex + non.european, subset(ukb.comb, date_1010 != 0), ties = "breslow")
summary(cox.tmp)

## Effects on outcomes seem to be attenuated in some disease constellations..
res.inter[ pval.inter < .05/nrow(res.inter) ]

##################################################
####        partial correlation matrix        ####
##################################################

## import partial correlation networks across all, women, and men
pcor.network        <- lapply(c("Both", "Male", "Female"), function(x){
  ## import the relevant network
  tmp     <- fread(paste0("../output/partial.correlation.", x,".UKB.phecodes.COVID19.20230904.txt"))
  ## add population
  tmp$pop <- x
  return(tmp)
})
## combine
pcor.network        <- do.call(rbind, pcor.network)
## edit names
names(pcor.network) <- c("id.1", "id.2", "estimate", "pval", "cor", "pop")
pcor.network        <- as.data.table(pcor.network)
## reshape
pcor.network        <- dcast(pcor.network, id.1 + id.2 ~ pop, value.var = c("estimate", "pval", "cor"))

## add phecode labels and categories
pcor.network        <- merge(pcor.network, lab.phe[, c("id", "phenotype", "category", "sex")], by.x="id.1", by.y = "id")
pcor.network        <- merge(pcor.network, lab.phe[, c("id", "phenotype", "category", "sex")], by.x="id.2", by.y = "id", suffixes = c(".1", ".2"))

## create sex-aware columns
pcor.network[, sex.edge := paste(sex.1, sex.2, sep="-")]
## drop diseases that cannot be related to each other
pcor.network        <- pcor.network[ !(sex.edge %in% c("Male-Female", "Female-Male"))]

## create combined estimates
pcor.network[, estimate := ifelse(sex.1 == "Female" | sex.2 == "Female", estimate_Female, 
                                  ifelse(sex.1 == "Male" | sex.2 == "Male", estimate_Male, estimate_Both))]
## same for p-values
pcor.network[, pval := ifelse(sex.1 == "Female" | sex.2 == "Female", pval_Female, 
                                  ifelse(sex.1 == "Male" | sex.2 == "Male", pval_Male, pval_Both))]

## some parallel computing
registerDoMC(10)

## try to find optimal solution by testing for adherence to a power-law
pcor.law            <- lapply(seq(0,.2,.01), function(x){
  
  ## get the network
  tmp  <- subset(pcor.network, (pval < .05/nrow(pcor.network) & estimate > x))
  
  ## compute the network
  tmp  <- graph_from_data_frame(tmp, direct = F, vertices = lab.phe[ id %in% tmp$id.1 | id %in% tmp$id.2]) 
  
  ## node degrees
  ndeg <- degree(tmp)
  
  ## test whether they follow a power law (https://cran.r-project.org/web/packages/poweRlaw/vignettes/b_powerlaw_examples.pdf)
  plaw <- displ$new(ndeg)
  ## infer model parameters
  est  <- estimate_xmin(plaw)
  ## update object
  plaw$setXmin(est)
  ## derive p-value
  plaw <- bootstrap_p(plaw, threads = 10, no_of_sims = 500)
  
  ## report back
  return(data.frame(cutoff=x, n.edges=ecount(tmp), n.nodes=length(ndeg), pval.law=plaw$p))
  
})
## combine
pcor.law <- do.call(rbind, pcor.law)
## choose pcor>0.02 as a cut-off, based on a power law

## put into a graph
tmp                 <- subset(pcor.network, (pval < .05/nrow(pcor.network) & estimate > .02))
pcor.graph          <- graph_from_data_frame(tmp, direct = F, vertices = lab.phe[ id %in% tmp$id.1 | id %in% tmp$id.2]) 
vcount(pcor.graph) ## n = 1381
ecount(pcor.graph) ## n = 5212

## compute layout
l                   <- layout_with_fr(pcor.graph, grid = "nogrid")
l                   <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

pdf("../graphics/Phecode.partial.correlation.network.20230913.pdf", width = 6.3, height = 6.3)
par(mar=rep(1,4), tck=-.01, cex.axis=.4, cex.lab=.4, mgp=c(.6,0,0), bty="n", yaxs="i", xaxs="i", lwd=.1)

## plot the network
plot(pcor.graph, layout=l*.9, rescale=F, vertex.size=2, 
     vertex.color=V(pcor.graph)$cl, 
     vertex.label=NA, 
     edge.color=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network) & E(pcor.graph)$estimate > .02, colorspace::adjust_transparency("grey80", .8), NA),
     edge.width=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network), E(pcor.graph)$estimate*10, 0))
## frame
rect(-1.1,-1.1,1.1,1.1, border = "black", col=NA, lwd=.5, xpd=NA)

## add legend for phecodes
pm <- par("usr")
legend(pm[1]-(pm[2]-pm[1])*.01, pm[3]+(pm[4]-pm[3])*.05, lty=0, pch=22, pt.lwd=.1, bty="n", xpd=NA,
       pt.bg=cl.cat$cl, legend = stringr::str_to_sentence(cl.cat$category),
       cex=.5, ncol=6, pt.cex=1)

dev.off()

## --> graph statistics <-- ##

## create table
pcor.charact <- data.table(name=names(V(pcor.graph)),
                           ## degree
                           deg = degree(pcor.graph),
                           ## betweeness
                           bet = betweenness(pcor.graph),
                           ## centrality
                           central = closeness(pcor.graph, mode="all", weights=NA),
                           ## hub status
                           hub = hub_score(pcor.graph, weights=NA)$vector,
                           ## authority
                           authority = authority_score(pcor.graph, weights=NA)$vector)
## add label
pcor.charact <- merge(lab.phe, pcor.charact, by.x="id", by.y="name")

## --> perform community detection <-- ##

## perform detection
pcor.com                <- cluster_edge_betweenness(pcor.graph)

## add to the graph object
V(pcor.graph)$community <- pcor.com$membership

## add to node statistics
tmp                     <- membership(pcor.com)
tmp                     <- data.frame(id = names(tmp), community = as.vector(tmp))
pcor.charact            <- merge(pcor.charact, tmp)
## n = 31 communities

## define a new colour vector
cl                      <- rainbow(length(unique(pcor.charact$community)))

pdf("../graphics/Phecode.partial.correlation.network.community.20230913.pdf", width = 6.3, height = 6.3)
par(mar=rep(1,4), tck=-.01, cex.axis=.4, cex.lab=.4, mgp=c(.6,0,0), bty="n", yaxs="i", xaxs="i", lwd=.1)

## plot the network
plot(pcor.graph, layout=l*.9, rescale=F, vertex.size=2, 
     vertex.color=cl[V(pcor.graph)$community], 
     vertex.label=NA, 
     edge.color=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network) & E(pcor.graph)$estimate > .02, colorspace::adjust_transparency("grey80", .8), NA),
     edge.width=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network), E(pcor.graph)$estimate*10, 0))
## frame
rect(-1.1,-1.1,1.1,1.1, border = "black", col=NA, lwd=.5, xpd=NA)

## add legend for phecodes
pm <- par("usr")
legend(pm[1]-(pm[2]-pm[1])*.01, pm[3]+(pm[4]-pm[3])*.05, lty=0, pch=22, pt.lwd=.1, bty="n", xpd=NA,
       pt.bg=cl, legend = 1:length(cl),
       cex=.5, ncol=12, pt.cex=1)

dev.off()

#-------------------------------------#
##--         summary figure        --##
#-------------------------------------#

## create new layout incorporating community structure
edge.weights <- function(community, network, weight.within = 100, weight.between = 1) {
  bridges <- crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights)
}

## assign new weights to each edge
E(pcor.graph)$weight <- edge.weights(pcor.com, pcor.graph, weight.within = 5, weight.between = .5)
l.weights            <- layout_with_fr(pcor.graph, grid = "nogrid")
l.weights            <- norm_coords(l.weights, ymin=-1, ymax=1, xmin=-1, xmax=1)

## open device
pdf("../graphics/Summary.disease.network.20230928.A.pdf", width = 6.3, height = 2.08)
## plotting parameters
par(mar=rep(1.5,4), tck=-.01, cex.axis=.4, cex.lab=.4, mgp=c(.6,0,0), bty="n", yaxs="i", xaxs="i", lwd=.1, mfrow=c(1,3))

#---------------------------------#
##--  disease network/ ICD-10  --##
#---------------------------------#

## plot the network
plot(pcor.graph, layout=l.weights*.9, rescale=F, vertex.size=2, 
     vertex.color=V(pcor.graph)$cl, 
     vertex.frame.width=.1,
     vertex.label=NA, 
     edge.color=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network) & E(pcor.graph)$estimate > .02, colorspace::adjust_transparency("grey80", .8), NA),
     edge.width=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network), E(pcor.graph)$estimate*3, 0))

## add legend for phecodes
pm <- par("usr")
legend(pm[1]-(pm[2]-pm[1])*.01, pm[3]+(pm[4]-pm[3])*.05, lty=0, pch=22, pt.lwd=.1, bty="n", xpd=NA,
       pt.bg=cl.cat$cl, legend = c(stringr::str_to_sentence(cl.cat$category[-nrow(cl.cat)]), "Other"),
       cex=.5, ncol=8, pt.cex=1)

## add legend for dots and lines
legend(pm[1]+(pm[2]-pm[1])*.02, pm[4]-(pm[4]-pm[3])*.02, lty=c(1,0), pch=c(NA, 21), pt.lwd=.3, col="black", lwd=.5, pt.bg="grey50",
       legend = c("|partial correlation|", "disease"), bty="n", cex=.6)

## add header
mtext("Disease-disease network", cex=.5)
## letter
text(pm[1], pm[4], labels = "A", font = 2)

#------------------------------------#
##-- highlight cardio-respiratory --##
#------------------------------------#

## plot the network
plot(pcor.graph, layout=l.weights*.9, rescale=F, vertex.size=2, 
     vertex.color=ifelse(V(pcor.graph)$community == 1, V(pcor.graph)$cl, "white"), 
     vertex.frame.width=.1,
     vertex.label=NA, 
     edge.color=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network) & E(pcor.graph)$estimate > .02, colorspace::adjust_transparency("grey80", .8), NA),
     edge.width=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network), E(pcor.graph)$estimate*3, 0))

## add header
mtext("Cardio-respiratory community", cex=.5)
## letter
text(pm[1], pm[4], labels = "B", font = 2)

#------------------------------------#
##--  highlight endocrine-renal   --##
#------------------------------------#

## plot the network
plot(pcor.graph, layout=l.weights*.9, rescale=F, vertex.size=2, 
     vertex.color=ifelse(V(pcor.graph)$community == 13, V(pcor.graph)$cl, "white"), 
     vertex.frame.width=.1,
     vertex.label=NA, 
     edge.color=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network) & E(pcor.graph)$estimate > .02, colorspace::adjust_transparency("grey80", .8), NA),
     edge.width=ifelse(E(pcor.graph)$pval < .05/nrow(pcor.network), E(pcor.graph)$estimate*3, 0))

## add header
mtext("Endocrine-renal community", cex=.5)
## letter
text(pm[1], pm[4], labels = "C", font = 2)

dev.off()

## open device
# png("../graphics/Summary.disease.network.20230928.B.png", width = 16, height = 8, res = 900, units = "cm")
pdf("../graphics/Summary.disease.network.20230928.B.pdf", width = 6.3, height = 3.15)
## plotting parameters
par(mar=rep(1.5,4), tck=-.01, cex.axis=.4, cex.lab=.4, mgp=c(.6,0,0), bty="n", yaxs="i", xaxs="i", lwd=.1)
layout(matrix(1:5, 1, 5), widths = c(.4,.15,.15,.15,.15))

#---------------------------------#
##--        Hub status         --##
#---------------------------------#

## change plotting settings
par(mar=c(2.5,13,1,.5), mgp=c(1.4,0,0), cex.lab=.7, cex.axis=.7)

## get the 30 most important hubs
tmp <- pcor.charact[ order(-hub)]
tmp <- tmp[ 1:30 ]

## empty plot
plot(c(0,1), c(.5,nrow(tmp)+.5), type="n", xlab="Hub score\n", xaxt="n", ylab="", yaxt="n",
     ylim=rev(c(.5, nrow(tmp)+.5)))
## add x-axis
axis(1, lwd=.5)
## add scores
rect(0, 1:nrow(tmp)-.4, tmp$hub, 1:nrow(tmp)+.4, lwd=.1, col=tmp$cl)
## add labels
pm <- par("usr")
## shorten one
text(pm[1], 1:nrow(tmp), pos=2, cex=.5, offset=.1, labels = gsub("Malignant neoplasm of other and ill-defined sites within the digestive organs and peritoneum",
                                                                 "Malignant neoplasm within the digestive organs and peritoneum", tmp$phenotype), xpd=NA)
## letter
text(pm[1]-(pm[2]-pm[1])*2, pm[4], labels = "C", font = 2, xpd=NA)

#---------------------------------#
##--       COVID-19 stats      --##
#---------------------------------#

## change plotting settings
par(mar=c(2.5,.5,1,.5))

## loop over all outcomes
for(j in lab.covid$covid){
  
  ## get the stats needed
  # foo <- merge(tmp, res.w.primary[ outcome == j], by=names(tmp)[1:11])
  foo <- res.w.primary[ outcome == j & id %in% tmp$id]
  ## order by Hub status again
  foo <- foo[ order(-hub)]
  
  ## empty plot
  plot(c(-1,2.8), c(.5,nrow(tmp)+.5), type="n", xlab=paste0("Hazard ratio\n", lab.covid$label[which(lab.covid$covid == j)]), xaxt="n", ylab="", yaxt="n",
       ylim=rev(c(.5, nrow(tmp)+.5)))
  ## add x-axis
  axis(1, lwd=.5, at=log(c(.5,1,2,5,10)), labels = c(.5,1,2,5,10))
  ## add rectangles to divide
  pm <- par("usr")
  rect(pm[1], 1:nrow(tmp)-.5, pm[2], 1:nrow(tmp)+.5, border=NA, col = c("white", "grey80"))
  ## no effect
  abline(v=0, lwd=.5, lty=2, xpd=F)
  
  ## add confidence intervals
  arrows(foo$beta - 1.96*foo$se, 1:nrow(foo), foo$beta + 1.96*foo$se, 1:nrow(foo), length = 0, col="grey50", lwd=.5, xpd=F)
  ## add point estimates
  points(foo$beta, 1:nrow(foo), pch=22, cex=.8, lwd=.3, bg=ifelse(foo$pval < .05/nrow(res.w.primary), "black", "white"))
  
  ## rectangle surrounding
  rect(pm[1], pm[3], pm[2], pm[4], lwd=.5)
}

dev.off()

## --> enrichment in certain communities <-- ##

## convert to data table
pcor.charact <- as.data.table(pcor.charact)

## enrichment of certain categories
enric.com    <- as.data.frame(expand.grid(category=unique(pcor.charact$category), community=unique(pcor.charact$community)))
enric.com    <- lapply(1:nrow(enric.com), function(x){
  
  ## get outcome
  outc <- enric.com$category[x]
  ## category
  dcat <- enric.com$community[x]
  
  ## community and category
  d1  <- nrow(pcor.charact[ category == outc & community == dcat ])
  ## non-category but community
  d2  <- nrow(pcor.charact[ category != outc & community == dcat ])
  ## category but not community
  d3  <- nrow(pcor.charact[ category == outc & community != dcat ])
  ## neither
  d4  <- nrow(pcor.charact[ category != outc & community != dcat ])
  
  ## test for enrichment
  enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))

  ## return information needed
  return(data.frame(enric.com[x,], or=enr$estimate, pval=enr$p.value, n.com=nrow(pcor.charact[ community == dcat ]), n.cat=nrow(pcor.charact[ category == outc ])))
  
})
## combine
enric.com <- do.call(rbind, enric.com)
## makes sense, with most community enriched for at most one disease category, with few exceptions,

## write to file
write.table(enric.com, "Enrichment.disease.categories.network.communities.20230915.txt", sep="\t", row.names = F)

## store other results as well
write.table(subset(pcor.network, (pval < .05/nrow(pcor.network) & estimate > .02)), "UKB.phecode.partial.correlation.network.20230926.txt", sep="\t", row.names=F)
write.table(pcor.charact, "UKB.phecode.partial.correlation.network.nodes.20230926.txt", sep="\t", row.names=F)

#-------------------------------------#
##-- map findings onto the network --##
#-------------------------------------#

## open plotting device
png("../graphics/Summary.phecodes.COVID19.network.mapping.20230913.png", units = "cm", width = 16, height = 16, res = 600)

## plotting parameters
par(mar=rep(1,4), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, lwd=.05, tck=-.01, mfcol=c(2,2))

## loop over all outcomes
for(j in lab.covid$covid){
  
  ## define color vector
  col.assoc <- colorRampPalette(c("white", "red1"))(max(-log10(res.w.primary[ outcome == j]$pval), na.rm=T)+1)
  
  ## --> network w/ primary care data <-- ##
  
  ## plot the network
  plot(pcor.graph, layout=l*.9, rescale=F, vertex.size=2, 
       ## color by association strength with COVID-19 outcome
       vertex.color=col.assoc[sapply(V(pcor.graph)$name, function(x){
         ## get the value+
         if(length(which(res.w.primary$id == x & res.w.primary$outcome == j)) > 0){
           return(ceiling(-log10(res.w.primary$pval[which(res.w.primary$id == x & res.w.primary$outcome == j)])))
         }else{
           return(1)
         }
       })], 
       vertex.label=NA, 
       vertex.frame.width=.1,
       edge.color=colorspace::adjust_transparency("grey80", .8),
       edge.width=E(pcor.graph)$estimate)
  ## frame
  rect(-1.1,-1.1,1.1,1.1, border = "black", col=NA, lwd=.5, xpd=NA)
  ## add header
  mtext(lab.covid$label[which(lab.covid$covid == j)], cex=.5, line=.2)
  
  ## colour gradient
  pm <- par("usr")
  ## length
  ll <- seq(pm[1]+(pm[2]-pm[1])*.05, pm[1]+(pm[2]-pm[1])*.35, length.out = length(col.assoc))
  ## rectangle for the colours
  rect(ll-(ll[2]-ll[1])/2, pm[3]+(pm[4]-pm[3])*.05, ll+(ll[2]-ll[1])/2, pm[3]+(pm[4]-pm[3])*.1, border=NA, col=col.assoc)
  ## box
  rect(ll[1]-(ll[2]-ll[1])/2, pm[3]+(pm[4]-pm[3])*.05, ll[length(ll)]+(ll[2]-ll[1])/2, pm[3]+(pm[4]-pm[3])*.1, border="black", col=NA, lwd=.3)
  ## add header
  text(pm[1]+(pm[2]-pm[1])*.05, pm[3]+(pm[4]-pm[3])*.12, cex=.6, labels = expression(-log[10]("p-value")), pos=4,
       offset = .2)
  ## simple axis
  text(ll[round(c(1, c(.2, .4, .6, .8, 1)*length(col.assoc)))], pm[3]+(pm[4]-pm[3])*.04,
       labels=round(c(0, c(.2, .4, .6, .8, 1)*length(col.assoc))), pos=1, cex=.4, offset = .1)
  
}

## close device
dev.off()

## write network and layout to file
save(pcor.graph, file="UKB.disease.network.20231023.RData")
save(l.weights, file="layout.UKB.disease.network.20231023.RData")

#------------------------------------#
##-- relation to COVID-19 results --##
#------------------------------------#

## add community and other characteristics to results (careful, not all diseases included)
res.w.primary <- merge(res.w.primary, pcor.charact, all.x=T, by=names(pcor.charact)[1:11])
## same for wide
res.wide      <- merge(res.wide, pcor.charact, all.x = T, by=names(pcor.charact)[1:11])

## --> enrichment in certain communities <-- ##

## enrichment of certain categories
res.enric.com <- as.data.frame(expand.grid(covid=lab.covid$covid, community=unique(res.w.primary$community)))
res.enric.com <- lapply(1:nrow(res.enric.com), function(x){
  
  ## get outcome
  outc <- res.enric.com$covid[x]
  ## category
  dcat <- res.enric.com$community[x]
  
  ## sig diseases in category
  d1  <- nrow(res.w.primary[ outcome == outc & community == dcat & pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3])
  ## non-sig diseases in category
  d2  <- nrow(res.w.primary[ outcome == outc & community == dcat & pval > .05/nrow(res.w.primary) & p.resid.phecode > 1e-3])
  ## sig diseases not in category
  d3  <- nrow(res.w.primary[ outcome == outc & community != dcat & pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3])
  ## non-sig diseases not in category
  d4  <- nrow(res.w.primary[ outcome == outc & community != dcat & pval > .05/nrow(res.w.primary) & p.resid.phecode > 1e-3])
  
  ## test for enrichment
  enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
  
  ## return information needed
  return(data.frame(res.enric.com[x,], or=enr$estimate, pval=enr$p.value, n.com=nrow(res.w.primary[ outcome == outc & community == dcat ])))
  
})
## combine
res.enric.com <- do.call(rbind, res.enric.com)

## write to file
write.table(res.enric.com, "Enrichment.COVID19.network.communities.20230915.txt", sep="\t", row.names = F)

## how many communities enriched for severe COVID-19 associations
View(aggregate(pval ~ community, subset(res.enric.com, or > 0 & covid != "long.covid.status"), function(x) sum(x < .05/nrow(res.enric.com))))

## --> test whether COVID-19 associations correlate with hub status <-- ##

## loop over all outcomes
cor.hub <- lapply(1:nrow(lab.covid), function(x){
  ## test 
  ct <- cor.test(res.wide$hub, unlist(-log10(res.wide[, paste0("pval.", lab.covid$covid[x]), with=F])), m="s", u="p")
  ## adjust for case numbers
  lt <- summary(lm(paste0("-log10(pval.",lab.covid$covid[x], ") ~ hub + n.phecode.", lab.covid$covid[x]), res.wide))$coefficients
  ## return results
  return(data.frame(lab.covid[x,], estimate=ct$estimate, pval=ct$p.val, beta=lt[2,1], se=lt[2,2], pval.lin=lt[2,4]))
})
## combine
cor.hub <- do.call(rbind, cor.hub)

##################################################
####       create MM counts by community      ####
##################################################

## re-import the relevant data set
ukb.comb <- fread("UKB.phecodes.COVID19.20230904.txt")
ukb.comb <- as.data.frame(ukb.comb)

## compute new exposure variables
for(j in 1:31){
  ## create new variable
  ukb.comb[, paste0("mm_", j)] <- rowSums(ukb.comb[, pcor.charact$id[ which(pcor.charact$community == j)]])
}
## look at newly derived variables
summary(ukb.comb[, paste0("mm_", 1:31)])

#----------------------------------#
##--  run association testing   --##
#----------------------------------#

## do in parallel
registerDoMC(12)

## define grid to be used
res.mm <- merge(data.frame(exposure=grep("rand_|mm_", names(ukb.comb), value=T)), lab.covid)

## run testing
res.mm <- mclapply(1:nrow(res.mm), function(x){
  
  ## get exposure
  expo <- res.mm$exposure[x]
  ## get outcome
  outc <- res.mm$covid[x]
  ## get date column
  cdat <- res.mm$date[x]
  
  ## data set needed (depends on disease)
  adj  <- "age.covid + sex.bin + pop"
  tmp  <- as.data.table(ukb.comb[, c("age.covid", "sex.bin", "pop", "non.european", expo, outc, cdat)])
  
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
      tmp <- data.frame(outcome=outc, phecode=expo, beta=ff$coefficients[1,1], se=ff$coefficients[1,3], pval=ff$coefficients[1,5], nevent=ff$nevent, nall=sum(!is.na(tmp[, outc, with=F])),
                        p.resid.phecode=ff.p$table[expo, 3], p.resid.overall=ff.p$table["GLOBAL", 3])
      ## add estimates for each strata
      jj  <- grep("strata", rownames(ff.s$coefficients), value=T)
      for(j in 1:length(jj)){
        ## estimate
        tmp[, paste0("beta.", j, ".ti")] <- ff.s$coefficients[jj[j],1]
        ## standard error
        tmp[, paste0("se.", j, ".ti")] <- ff.s$coefficients[jj[j],2]
        ## estimate
        tmp[, paste0("pval.", j, ".ti")] <- ff.s$coefficients[jj[j],5]
      }
      ## return
      return(tmp)
    }
  }
}, mc.cores = 12) 
## combine into one
res.mm <- do.call(plyr::rbind.fill, res.mm)
## add z-score to enable sorting
res.mm <- as.data.table(res.mm)
res.mm[, zscore := abs(beta/se) ]
## HR column
res.mm[, hr_column := paste0(sprintf("%.2f", exp(beta)), " (", sprintf("%.2f", exp(beta - 1.96*se)),";", sprintf("%.2f", exp(beta + 1.96*se)),")")]

##################################################
####      convergence genetic correlation     ####
##################################################

## import genetic correlation
res.cor       <- fread("<path to file>/covid_release7_vs_allphecodes.LF.20221214")
## Long Covid (strict case and broad control definition)
tmp           <- fread("<path to file>/all_results_20231909.txt")
## keep only one endpoint (according to preprint)
tmp           <- tmp[ p1 == "st_cas_vs_bd_con.sumstats.gz"]
## edit some names
tmp[, p1 := gsub("\\.sumstats\\.gz", "", p1)]
tmp[, p2 := gsub("\\.sumstats\\.gz", "", p2)]

## use same population as for the MR
res.cor       <- rbind(res.cor[ p1 %in% c("A2_ALL_exUKBB", "B2_ALL_exUKBB")], tmp)
## delete missing values
res.cor       <- res.cor[ !is.na(rg) ]
## how many phecodes
length(unique(res.cor$p2))
## n = 1200 (n=3546 in total)

## create bridging file
tmp.brd        <- data.frame(covid=c("hosp.covid19", "severe_resp", "death.covid19", "long.covid.status"), 
                             hgi=c("B2_ALL_exUKBB", "A2_ALL_exUKBB", "A2_ALL_exUKBB", "st_cas_vs_bd_con"))
## add to genetic correlation results (careful, duplicates A2)
res.cor        <- merge(res.cor, tmp.brd, by.x="p1", by.y="hgi", allow.cartesian = T)
## add survival estimates
res.cor        <- merge(res.cor, res.w.primary, by.x=c("covid", "p2"), by.y=c("outcome", "id"))
## looses some with no observational estimates

## compute convergence; use p-value cut-off for number of associations in obs results
res.cor[, convergence := ifelse(pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3 & p < .05/1128 & sign(rg) == sign(beta), "yes", "no")]
table(res.cor$convergence) ## 57 unique endpoints
#   no  yes 
# 3940   75 
plot(rg ~ beta, res.cor[ covid == "hosp.covid19" & pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3])

## export results
fwrite(res.cor[pval < .05/nrow(res.w.primary) & p.resid.phecode > 1e-3], "Results.survival.genetic.correlation.txt", sep="\t", row.names=F)

#----------------------------------#
##--           Figure           --##
#----------------------------------#

## open device
pdf("../graphics/Convergence.survival.genetic.correlation.20230918.pdf", width = 6.3, height = 6.3)
## graphical parameters
par(mar=c(2.5,13,.5,.5), cex.axis=.6, cex.lab=.6, mgp=c(1,0,0), tck=-.01, xaxs="i", yaxs="i")
## define layout
layout(matrix(1:5, 1, 5), widths = c(.4,rep(.15,4)))

## tmp data for plotting
tmp    <- lab.cox[ id %in% res.cor[ convergence == "yes"]$p2 ]
## order 
tmp    <- tmp[ order(srt) ]
## new sorting
tmp[, plt.srt := 1:nrow(tmp)]
## res.cor to be plotted
tmp    <- merge(res.cor, tmp[, c("id", "plt.srt")], by.x="p2", by.y="id", all.y=T)

## colouring scheme
cl.tmp <- do.call(data.frame, aggregate(plt.srt ~ category, unique(tmp[, c("plt.srt", "category")]), function(x) c(min(x), mean(x), max(x))))
## add colour
cl.tmp <- merge(cl.tmp, cl.cat)
cl.tmp <- cl.tmp[ order(cl.tmp$plt.srt.1), ]

#-----------------------------#
##--     Hazard ratios     --##
#-----------------------------#

## loop through severe endpoints
for(j in lab.covid$covid[1:3]){
  
  ## adapt plotting parameters
  if(j != "hosp.covid19") par(mar=c(2.5,.5,.5,.5))
  
  ## subset of outcome to be plotted
  foo <- tmp[covid == j]
  
  ## empty plot
  plot(c(-1,2.8), c(.5,nrow(foo)+.5), type="n", xlab=paste0("Hazard ratio\n", lab.covid$label[which(lab.covid$covid == j)]), xaxt="n", ylab="", yaxt="n",
       ylim=rev(c(.5, nrow(foo)+.5)))
  ## add x-axis
  axis(1, lwd=.5, at=log(c(.5,1,2,5,10)), labels = c(.5,1,2,5,10))
  ## add rectangles to divide
  pm <- par("usr")
  # rect(pm[1], 1:nrow(foo)-.5, pm[2], 1:nrow(foo)+.5, border=NA, col = c("white", "grey80"))
  rect(pm[1], cl.tmp$plt.srt.1-.5, pm[2], cl.tmp$plt.srt.3+.5, border=NA, col=adjust_transparency(cl.tmp$cl, .5))
  ## no effect
  abline(v=0, lwd=.5, lty=2, xpd=F)
  
  ## add confidence intervals (se.y refers to estimate from survival models) 
  arrows(foo$beta - 1.96*foo$se.y, foo$plt.srt, foo$beta + 1.96*foo$se.y, foo$plt.srt, length = 0, col="grey50", lwd=.5, xpd=F)
  ## add point estimates
  points(foo$beta, foo$plt.srt, pch=22, cex=.8, lwd=.3, bg=ifelse(foo$pval < .05/nrow(res.w.primary), "black", "white"))
  
  ## rectangle surrounding
  rect(pm[1], pm[3], pm[2], pm[4], lwd=.5)
  
  ## add labels
  if( j == "hosp.covid19") text(pm[1], foo$plt.srt, labels = foo$phenotype, cex=.6, xpd=NA, pos=2, offset=.1)
}

#-----------------------------#
##--  Genetic correlation  --##
#-----------------------------#

## loop through severe endpoints
for(j in c("B2", "A2")){
  
  ## subset of outcome to be plotted
  foo <- tmp[covid %in% c("hosp.covid19", "death.covid19") & p1 == paste0(j, "_ALL_exUKBB")]
  
  ## empty plot
  plot(c(-.1,.8), c(.5,nrow(foo)+.5), type="n", xlab=paste0("Genetic correlation\n", ifelse(j == "B2", "Hospitalisation", "Critical illness")), xaxt="n", ylab="", yaxt="n",
       ylim=rev(c(.5, nrow(foo)+.5)))
  ## add x-axis
  axis(1, lwd=.5, at=c(0,.25,.5,.75))
  ## add rectangles to divide
  pm <- par("usr")
  # rect(pm[1], 1:nrow(foo)-.5, pm[2], 1:nrow(foo)+.5, border=NA, col = c("white", "grey80"))
  rect(pm[1], cl.tmp$plt.srt.1-.5, pm[2], cl.tmp$plt.srt.3+.5, border=NA, col=adjust_transparency(cl.tmp$cl, .5))
  ## no effect
  abline(v=0, lwd=.5, lty=2, xpd=F)
  
  ## add confidence intervals (se.x refers to estimate from genetic correlation) 
  arrows(foo$rg - 1.96*foo$se.x, foo$plt.srt, foo$rg + 1.96*foo$se.x, foo$plt.srt, length = 0, col="grey50", lwd=.5, xpd=F)
  ## add point estimates
  points(foo$rg, foo$plt.srt, pch=22, cex=.8, lwd=.3, bg=ifelse(foo$p < .05/1128, "black", "white"))
  
  ## rectangle surrounding
  rect(pm[1], pm[3], pm[2], pm[4], lwd=.5)
  
}
## close device
dev.off()

#-----------------------------#
##--   import MR results   --##
#-----------------------------#

## import mr results
res.mr <- fread("<path to file>/Results.MR.phecode.to.COVID19.20221213.txt")
## create name to ease merging with correlation results
res.mr[, p1 := paste0(pheno, "_ALL_exUKBB")]
## subset to results with evidence of convergence
res.mr <- merge(res.mr, res.cor, by.x = c("id", "phecode", "phenotype", "category", "p1"), by.y = c("p2", "phecode", "phenotype", "category", "p1"))
## subset to convergence
res.mr <- res.mr[ convergence == "yes" ]
## how many positive findings
res.mr[ P.value.presso.combined < .05 ]
## create MR column
res.mr[, beta_column := paste0(sprintf("%.2f", ifelse(!is.na(P.value.presso.wo), exp(Estimate.presso.wo), exp(Estimate.presso))), " (",
                               sprintf("%.2f", ifelse(!is.na(P.value.presso.wo), exp(Estimate.presso.wo - 1.96*Std.Error.presso.wo), exp(Estimate.presso - 1.96*Std.Error.presso))), "; ",
                               sprintf("%.2f", ifelse(!is.na(P.value.presso.wo), exp(Estimate.presso.wo + 1.96*Std.Error.presso.wo), exp(Estimate.presso + 1.96*Std.Error.presso))),")")]

## export for tables
fwrite(res.mr, "Results.convergence.genetic.correlation.mr.covid19.20230918.txt", sep="\t", row.names = F)

##################################################
####  Supplemental Tables for the manuscript  ####
##################################################

## import phecode label from genetic analysis
lab.gen <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/14_phecode_GWAS/07_summary_discovery/data/Results.LDSC.phecode.GWAS.with.mapping.previous.GWAS.curated.20221215.txt")
## to have all codes in one place
tmp     <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/10_phecode_GWAS_mapping/data/Phecodes.GWAS.catalog.mapping.20211115.txt")

## create new file
st1     <- merge(lab.phe[, c("id", "category", "phecode", "phenotype", "sex", "cases.w.primary", "cases.wo.primary")],
                 lab.gen[, c("id", "cases", "controls", "primary_care", "wo.primary_care")], all.x=T)
st1     <- merge(st1, tmp[, c("id", "icd10", "read2", "ctv3", "snomed")])

## write to file
write.table(st1, "Supplementary.Table.1.case.counts.20231002.txt", sep="\t", row.names=F)

#------------------------------#
##-- cohort characteristics --##
#------------------------------#

## gender
table(ukb.comb$sex)/nrow(ukb.comb)
## age
mean(ukb.comb$age.covid); sd(ukb.comb$age.covid)
## ancestry
table(ukb.comb$non.european)/nrow(ukb.comb)

########################################################################################
########################################################################################
############          REVISION Communications Medicine 02/01/2024           ############
########################################################################################
########################################################################################


##################################################
####          multivariable adjusted          ####
##################################################

## obtain additional baseline confounders, including smoking, alcohol, BMI, SES
cl.select       <- c("f.eid", "f.20116.0.0", "f.1558.0.0", "f.189.0.0")
## names to assign
cl.names        <- c("f.eid", "smoking", "alcohol", "townsend")

## import data from the main release
ukb.dat         <- read_parquet("<path to file>/ukb45268.parquet",
                                col_select = cl.select)
## change names
names(ukb.dat)  <- cl.names

## define factor label for smoking
ukb.dat$smoking <- as.character(ukb.dat$smoking)
ukb.dat$smoking <- ifelse(ukb.dat$smoking == "Prefer not to answer", "Previous", ukb.dat$smoking)
ukb.dat$smoking <- factor(ukb.dat$smoking, levels = c("Never", "Previous", "Current"))

## define factor level for alcohol
ukb.dat$alcohol <- as.character(ukb.dat$alcohol)
ukb.dat$alcohol <- ifelse(ukb.dat$alcohol == "Prefer not to answer", "Three or four times a week", ukb.dat$alcohol)
ukb.dat$alcohol <- factor(ukb.dat$alcohol, levels = c("Never", "Special occasions only", "One to three times a month", "Once or twice a week", "Three or four times a week",
                                                      "Daily or almost daily"))
## convert to data table
ukb.dat         <- as.data.table(ukb.dat)

## health care utilization (number HES entries five years before pandemic)
hesin_diag      <- fread("<path to file>/hesin_diag.txt", sep="\t")
hesin           <- fread("<path to file>/hesin.txt", sep="\t")
## subset to 5 years before the pandemic
hesin           <- hesin[ as.Date(epistart, format="%d/%m/%Y") >= as.Date("01/01/2015", format="%d/%m/%Y") & as.Date(epistart, format="%d/%m/%Y") <= as.Date("01/01/2020", format="%d/%m/%Y")]
## create variables (how many occasions, and how many days spend in hospital)
hesin           <- hesin[, .(count = .N, days_hospital = sum(as.Date(epiend, format="%d/%m/%Y") - as.Date(epistart, format="%d/%m/%Y"))), by = eid]

## same for primary care utilization (use TPP and EMIS only)
ukb.tpp         <- fread("<path to file>/covid19_tpp_gp_clinical.txt")
## reduce to 5 years before pandemic and delete possible coding errors
ukb.tpp         <- ukb.tpp[ as.Date(event_dt, format="%d/%m/%Y") >= as.Date("01/01/2015", format="%d/%m/%Y") & as.Date(event_dt, format="%d/%m/%Y") <= as.Date("01/01/2020", format="%d/%m/%Y") ]
## create variables (how many occasions seen the GP)
ukb.tpp         <- ukb.tpp[, .(count = length(unique(event_dt))), by = eid]
## import EMIS data
ukb.emis        <- fread("<path to file>/covid19_emis_gp_clinical.txt")
## reduce to 5 years before pandemic and delete possible coding errors
ukb.emis        <- ukb.emis[ as.Date(event_dt, format="%d/%m/%Y") >= as.Date("01/01/2015", format="%d/%m/%Y") & as.Date(event_dt, format="%d/%m/%Y") <= as.Date("01/01/2020", format="%d/%m/%Y") ]
## create variables (how many occasions seen the GP)
ukb.emis        <- ukb.emis[, .(count = length(unique(event_dt))), by = eid]
## merge both (some participants overlap)
ukb.primary     <- merge(ukb.emis, ukb.tpp, by="eid", suffixes = c(".emis", ".tpp"), all=T)
## create novel count variable
ukb.primary[, count := ifelse(!is.na(count.emis), count.emis, 0) + ifelse(!is.na(count.tpp), count.tpp, 0) ]

## add to combined data set
ukb.comb        <- merge(ukb.comb, hesin, all.x=T, by.x="f.eid", by.y="eid")
ukb.comb        <- merge(ukb.comb, ukb.primary[, c("eid", "count"), with=F], all.x=T, by.x="f.eid", by.y="eid", suffixes = c(".hes", ".primary"))
ukb.comb        <- as.data.table(ukb.comb)
## recode NAs
ukb.comb[, count.hes := ifelse(is.na(count.hes), 0, count.hes)]
ukb.comb[, days_hospital := ifelse(is.na(days_hospital), 0, days_hospital)]
ukb.comb[, count.primary := ifelse(is.na(count.primary), 0, count.primary)]
## add other covariates
ukb.comb        <- merge(ukb.comb, ukb.dat, by="f.eid", all.x = T)

## additional diseases at baseline; including any two of alcohol abuse, anemia, arrhythmia, asthma, cancer, 
##                                  chronic kidney disease, chronic pulmonary disorders, cirrhosis, coagulopathy, 
##                                  congestive heart failure, chronic obstructive pulmonary disease, coronary artery disease, dementia, 
##                                  diabetes type 1, diabetes type 2, end-stage renal disease on dialysis, hemiplegia, HIV, hypertension, 
##                                  hypertension and type 1 or 2 diabetes diagnosis, inflammatory bowel disorder, lupus or systemic lupus erythematosus, 
##                                  mental health disorders, multiple sclerosis, Parkinson’s disease, peripheral vascular disorders, pregnant, 
##                                  pulmonary circulation disorder, rheumatoid arthritis, seizure/epilepsy, severe obesity (BMI > = 40 kg/m2), weight loss,
##                                  Down’s syndrome, other substance abuse, cystic fibrosis, autism, sickle cell, corticosteroid drug prescriptions, immunosuppressant drug prescriptions.
phecodes.conf   <- c(317, 317.1, 317.11, 280.1, 280.2, 281.11, 281.12, 281.13, 281.9, 282.9, 282.9, 283, 283.1, 284, 285, 285.1, 285.2, 285.22, 427.5,
                     427.2, 427.21, 427.22, 427.8, 495, 495.2, lab.phe$phecode[which(lab.phe$category == "neoplasms")], 401.22, 585.3, 585.4, 496, 496.2, 
                     496.21, 571.51, 571.5, 571.8, 286.7, 286.5, 286.11, 286.12, 286.13, 286.6, 428.1, 411.1, 411.4, 290, 290.1, 290.12, 290.16,
                     290.11, 250.1, 250.11, 250.12, 250.13, 250.13, 250.2, 250.21, 250.22, 250.22, 250.23, 250.24, 585.32, 585.1, 585.2, 585.31,
                     342, 71, 71.1, 401.1, 402, 415.21, 555.2, 555.21, 695.41, 695.42, 292.4, 295.1, 295.2, 295.3, 296, 296.1, 296.22, 297.2, 300,
                     300.1, 300.11, 300.12, 300.13, 300.3, 300.9, 301, 301.1, 301.2, 303.1, 303.3, 303.4, 305.2, 305.21, 306, 312, 312.3, 313.3,
                     335, 332, 443.8, 443.9, 443.1, 395.4, 415, 415.11, 415.2, 415.21, 714, 714.1, 714.2, 345, 345.1, 345.11, 345.12,
                     316, 499)
## make unique
phecodes.conf   <- unique(phecodes.conf)

## create variable whether at least two entries
ukb.comb[, morb_conf := apply(ukb.comb[, paste0("date_", phecodes.conf), with=F], 1, function(x) sum(x))]
table(ukb.comb$morb_conf)
## categorize
ukb.comb[, morb_conf := ifelse(morb_conf > 1 , 1, 0)]
table(ukb.comb$morb_conf)
#      0      1 
# 151546 287371 

## delete what is no longer needed
rm(ukb.emis); rm(ukb.tpp); rm(hesin); rm(hesin_diag); rm(ukb.primary); gc(reset=T)

## write to file for association testing
fwrite(ukb.comb, "UKB.phecodes.COVID19.20240102.txt", sep="\t", na = NA, row.names = F)

#-----------------------------------#
##--   augmented baseline model  --##
#-----------------------------------#

## Hospitalisation
cox.hospital.rev    <- coxph(Surv(fol.hosp, hosp.covid19) ~ age.covid + sex + pop + smoking + alcohol + bmi + count.hes + count.primary + days_hospital + morb_conf + townsend, 
                             ukb.comb, ties = "breslow")
summary(cox.hospital.rev)
cox.zph(cox.hospital.rev)$table[,3] ## bmi, count primary, and townsend with strong violation
## stratify by time interval
cox.hospital.rev    <- survSplit(Surv(fol.covid, hosp.covid19) ~ age.covid + sex + pop + smoking + alcohol + bmi + count.hes + count.primary + days_hospital + morb_conf + townsend,
                                 ukb.comb, cut=c(.5,1,1.5,2), episode = "tgroup")
## fit the model (need to redefine names of variables)
cox.hospital.rev    <- coxph(Surv(tstart, fol.covid, hosp.covid19) ~ age.covid:strata(tgroup) + sex + pop + smoking + alcohol + bmi + count.hes + count.primary + days_hospital + morb_conf + townsend,
                             cox.hospital.rev , ties = "breslow")
summary(cox.hospital.rev)
## associations as expected, with age, male sex, and non-European ancestry being strongly associated

## Severe respiratory failure
cox.severe.resp.rev  <- coxph(Surv(fol.covid, severe_resp) ~ age.covid + sex + pop + smoking + alcohol + bmi + count.hes + count.primary + days_hospital + morb_conf + townsend, 
                          ukb.comb, ties = "breslow")
summary(cox.severe.resp.rev)
cox.zph(cox.severe.resp.rev)$table[,3]

## COVID-19 death
cox.covid19.death.rev <- coxph(Surv(fol.death, death.covid19) ~ age.covid + sex + non.european + smoking + alcohol + bmi + count.hes + count.primary + days_hospital + morb_conf + townsend, 
                           ukb.comb, ties = "breslow")
summary(cox.covid19.death.rev)
cox.zph(cox.covid19.death.rev)$table[,3]

## Long-COVID
cox.long.covid.rev    <- coxph(Surv(fol.long.covid, long.covid.status) ~ age.covid + sex + pop + smoking + alcohol + bmi + count.hes + count.primary + days_hospital + morb_conf + townsend, 
                               ukb.comb, ties = "breslow")
summary(cox.long.covid.rev)
cox.zph(cox.long.covid.rev)
## looks okish

#-----------------------------------#
##--        import results       --##
#-----------------------------------#

## import results with additional adjustment
res.w.primary.adj <- fread("../output/Results.adjusted.UKB.phecodes.COVID19.20240102.txt")

## combine with crude results (some models did no longer converge)
res.w.primary.adj <- merge(res.w.primary, res.w.primary.adj, by.x=c("id", "outcome"), by.y=c("phecode", "outcome"), suffixes = c(".crude", ".adjusted"), all.x=T)

## overall attenuation
plot(-log10(pval.adjusted) ~ I(-log10(pval.crude)), res.w.primary.adj)

## how many findings are still passing statistical significance (careful, this includes 'newcomers')
nrow(res.w.primary.adj[ pval.adjusted < .05/nrow(res.w.primary.adj)])
## n = 668

## how many from those that were previously significant
nrow(res.w.primary.adj[  pval.crude < .05/nrow(res.w.primary) & p.resid.phecode.crude > 1e-3 & pval.adjusted < .05/1128])
## n = 718

## add coloumn to ease reporting
res.w.primary.adj[, hr_column.crude := paste0(sprintf("%.2f", exp(beta.crude)), " (", sprintf("%.2f", exp(beta.crude - 1.96*se.crude)),";", sprintf("%.2f", exp(beta.crude + 1.96*se.crude)),")")]
res.w.primary.adj[, hr_column.adjusted := paste0(sprintf("%.2f", exp(beta.adjusted)), " (", sprintf("%.2f", exp(beta.adjusted - 1.96*se.adjusted)),";", sprintf("%.2f", exp(beta.adjusted + 1.96*se.adjusted)),")")]

## write to file
write.table(res.w.primary.adj[  pval.crude < .05/nrow(res.w.primary) & p.resid.phecode.crude > 1e-3 & pval.adjusted < .05/1128], 
            "Results.additional.adjustment.UKB.phecodes.20240103.txt", sep="\t", row.names = F)

##################################################
####            interaction testing           ####
##################################################

## create new age binary
ukb.comb[, age.cat := ifelse(age.covid >= 65, 1, 0)]
table(ukb.comb$age.cat)
#      0      1 
# 157252 281665 

## create new townsend binary
ukb.comb[, townsend.cat := ifelse(townsend >= median(ukb.comb$townsend, na.rm = T), 1, 0)]
table(ukb.comb$townsend.cat, useNA = "always")
#        0      1   <NA> 
#   219178 219194    545 

## write to file for association testing
fwrite(ukb.comb, "UKB.phecodes.COVID19.20240102.txt", sep="\t", na = NA, row.names = F)

#-----------------------------------#
##--        import results       --##
#-----------------------------------#

## import results for possible effect modification
res.inter.rev <- lapply(c("sex", "non.european", "age.cat", "townsend.cat"), function(x){
  
  ## import results
  if(x %in% c("age.cat", "townsend.cat")){
    tmp <- fread(paste0("../output/Results.", x,".UKB.phecodes.COVID19.20240102.txt"))
  }else{
    tmp <- fread(paste0("../output/Results.", x,".UKB.phecodes.COVID19.20230904.txt"))
  }
  
  ## add interaction type
  tmp[, term := x]
  ## return
  return(tmp)
  
})
## combine
res.inter.rev <- rbindlist(res.inter.rev)
## add label
res.inter.rev <- merge(lab.phe, res.inter.rev, by.x = "id", by.y = "phecode")
res.inter.rev <- as.data.table(res.inter.rev)

## significant findings
res.inter.rev[ pval.inter < .05/nrow(res.inter.rev) ]
## n = 16

## look into one strong example
table(ukb.comb$date_275.3, ukb.comb$age.cat)
## sanity check 
cox.tmp <- coxph(Surv(fol.death, death.covid19) ~  date_275.3*age.cat + sex + non.european, ukb.comb, ties = "breslow")
summary(cox.tmp)
cox.zph(cox.tmp)$table[,3]
cox.tmp <- coxph(Surv(fol.death, death.covid19) ~  date_275.3 + sex + non.european, subset(ukb.comb, age.cat == 0), ties = "breslow")
summary(cox.tmp)
cox.zph(cox.tmp)$table[,3]
cox.tmp <- coxph(Surv(fol.death, death.covid19) ~  date_275.3 + sex + non.european, subset(ukb.comb, age.cat == 1), ties = "breslow")
summary(cox.tmp)
cox.zph(cox.tmp)$table[,3]

## create HR columns
res.inter.rev[, hr_column.1 := paste0(sprintf("%.2f", exp(beta.1)), " (", sprintf("%.2f", exp(beta.1 - 1.96*se.1)),";", sprintf("%.2f", exp(beta.1 + 1.96*se.1)),")")]
res.inter.rev[, hr_column.2 := paste0(sprintf("%.2f", exp(beta.2)), " (", sprintf("%.2f", exp(beta.2 - 1.96*se.2)),";", sprintf("%.2f", exp(beta.2 + 1.96*se.2)),")")]

## write to file
write.table(res.inter.rev[ pval.inter < .05/nrow(res.inter.rev)], "Sig.Results.phecode.COVID19.sex.euro.age.townsend.interaction.20240103.txt", sep="\t", row.names = F)

##################################################
####        look into competing risk          ####
##################################################

## create new variables that define 'death w/o COVID-19' as a competing risk
ukb.comb[, death.covid19.comprsk := ifelse( death.covid19 == 0 & !is.na(date_of_death), 2, death.covid19)]
table(ukb.comb$death.covid19.comprsk)
#      0      1      2 
# 424142   1546  13229 

## hospitalization
ukb.comb[, hosp.covid19.comprsk := ifelse( hosp.covid19 == 0 & !is.na(date_of_death), 2, hosp.covid19)]
table(ukb.comb$hosp.covid19.comprsk)
#      0      1      2 
# 418746   7507  12664 

## severe respiratory failure
ukb.comb[, severe_resp.comprsk := ifelse( severe_resp == 0 & !is.na(date_of_death), 2, severe_resp)]
table(ukb.comb$severe_resp.comprsk)
#      0      1      2 
# 423813    662  14442 

## Long COVID
ukb.comb[, long.covid.status.comprsk := ifelse( long.covid.status == 0 & !is.na(date_of_death), 2, long.covid.status)]
table(ukb.comb$long.covid.status.comprsk)
#      0      1      2 
# 423681    470  14765 

## write to file for data analysis
fwrite(ukb.comb, "UKB.phecodes.COVID19.20240102.txt", sep="\t", na = NA, row.names = F)

#-----------------------------------#
##--      import everything      --##
#-----------------------------------#

## import results from competing risk models
res.comprsk       <- fread("../output/Results.competing.UKB.phecodes.COVID19.20240102.txt")
## adapt names to merge with other results
res.comprsk[, outcome := gsub("\\.comprsk", "", outcome)]
 
## add to adjusted models to simplify reporting
res.w.primary.adj <- merge(res.w.primary.adj, res.comprsk, by.x=c("id", "outcome"), by.y=c("phecode", "outcome"), all.x=T)

## how many remained significant in comp risk model
nrow(res.w.primary.adj[pval.crude < .05/nrow(res.w.primary) & p.resid.phecode.crude > 1e-3 & pval < .05/1128])
## n = 1126

## simple plotting
plot(-log10(pval) ~ I(-log10(pval.crude)), res.w.primary.adj)
res.w.primary.adj[pval.crude < .05/nrow(res.w.primary) & p.resid.phecode.crude > 1e-3 & pval > .05/1128]
## what drops looks highly plausible

## add coloumn to ease reporting
res.w.primary.adj[, hr_column.comprsk := paste0(sprintf("%.2f", exp(beta)), " (", sprintf("%.2f", exp(beta - 1.96*se)),";", sprintf("%.2f", exp(beta + 1.96*se)),")")]

## add table with better naming of outcomes
res.w.primary.adj <- merge(res.w.primary.adj, lab.covid[, c("covid", "label")], by.x="outcome", by.y="covid")

## write to file (last columns refer to competing risk estimates)
write.table(res.w.primary.adj[  pval.crude < .05/nrow(res.w.primary) & p.resid.phecode.crude > 1e-3 ], 
            "Results.additional.adjustment.competing.risk.UKB.phecodes.20240103.txt", sep="\t", row.names = F)