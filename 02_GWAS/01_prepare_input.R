######################################################
#### GWAS on phecodes from UKBB EHR-mapping       ####
#### Maik Pietzner                                ####
######################################################

rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)
load(".RData")

## faster data handling
require(data.table)

##################################################
####      import covariate data variables     ####
##################################################

## import main UKBB release
load("<path to file>/ukb45268.RData")

## import data dictionary 
lab.main        <- read.table("<path to file>/Data.dictionary.UKBB.main.dataset.45268.txt",
                              sep="\t", header=T)

## retain only information needed
ukbb.dat        <- bd[, c("f.eid", "f.21022.0.0", "f.31.0.0", "f.54.0.0", "f.22000.0.0", paste0("f.22009.0.", 1:10))]
## rename
names(ukbb.dat) <- c("f.eid", "age", "sex", "centre", "batch", paste0("pc", 1:10))
## convert to data table
ukbb.dat        <- as.data.table(ukbb.dat)

## transform batch to binary
ukbb.dat$batch  <- ifelse(ukbb.dat$batch < 0, 1,  0)
table(ukbb.dat$batch)
rm(bd); gc()

## create dummy participant list to allow for inclusion in phecode data
tmp             <- ukbb.dat[, "f.eid", with = F]

## drop participants with missing values
ukbb.dat        <- na.omit(ukbb.dat)

## make sex binary
ukbb.dat$sex    <- ifelse(ukbb.dat$sex == "Female", 1, 0)

## write to file
fwrite(ukbb.dat, "covariates.txt", sep = "\t", row.names=F, quote=F)

## create second covariate file excluding sex from the adjustment set
fwrite(ukbb.dat[ sex == 1, c("FID", "IID", "age", "centre", "batch", paste0("pc", 1:10)), with=F], "covariates.female.txt", sep = "\t", row.names=F, quote=F)
fwrite(ukbb.dat[ sex == 0, c("FID", "IID", "age", "centre", "batch", paste0("pc", 1:10)), with=F], "covariates.male.txt", sep = "\t", row.names=F, quote=F)

##################################################
####        import phecode definitions        ####
##################################################

## import phecode data
require(data.table)
ukbb.phe     <- fread("<path to file>", sep="\t", header = T)

## convert to wide format
ukbb.phe     <- dcast(ukbb.phe, f.eid ~ phecode, value.var = c("date", "resource"))

## import label
lab.phe      <- read.table("ç", sep="\t", header=T)
## add identifier to match with data set
lab.phe$id   <- paste0("date_", lab.phe$phecode) 
## drop all with less than 100 cases
lab.phe      <- subset(lab.phe, n.charite >= 100) 

## add sex-specific coding
lab.tmp      <- read.csv("ç/phecode_definitions1.2.csv")
## create variable to map to UKBB data set
lab.tmp$id   <- paste0("date_", lab.tmp$phecode) 
## add sex-specific information to the data
lab.phe      <- merge(lab.phe, lab.tmp[, c("id", "sex")], by="id")
## replace missing ones
lab.phe$sex[lab.phe$sex == ""] <- "Both"

## add to covariates
ukbb.phe     <- merge(tmp, ukbb.phe, all.x=T, by="f.eid")

## convert dates to binary variables
ukbb.phe               <- as.data.frame(ukbb.phe)
ukbb.phe[, lab.phe$id] <- apply(ukbb.phe[, lab.phe$id], 2, function(x){
  x <- ifelse(!is.na(x), 1, 0)
  return(x)
})

#--------------------------------------#
##-- store batches of 100 diseases  --## 
#--------------------------------------#

## distinguish by sex-combined and sex-specific analysis
lab.phe            <- lab.phe[order(lab.phe$sex, lab.phe$phecode),]
## create indicator in the label data
lab.phe$batch      <- NA
## now sort again
lab.phe            <- lab.phe[order(lab.phe$sex, lab.phe$batch, lab.phe$phecode),]
## create additional batches
lab.phe$batch      <- c(rep(1:12, each=100), rep(13, 109), rep(14, 113), rep(15,26))

## store phenotypic data sets
for(j in 1:15){
  
  if(j %in% 1:13){
    ## store the subset of phenotypes in the current batch
    fwrite(ukbb.phe[, c("FID", "IID", lab.phe$id[which(lab.phe$batch == j)]), with = F],  paste0("phenotype.",j,".txt"),sep = "\t", row.names=F, quote=F, na = "NA")
  }else if(j == 14){
    ## only women
    fwrite(ukbb.phe[IID %in% subset(ukbb.dat, sex == 1)$IID, c("FID", "IID", lab.phe$id[which(lab.phe$batch == j)]), with = F],  paste0("phenotype.",j,".txt"),sep = "\t", row.names=F, quote=F, na = "NA")
  }else{
    ## only men
    fwrite(ukbb.phe[IID %in% subset(ukbb.dat, sex == 0)$IID, c("FID", "IID", lab.phe$id[which(lab.phe$batch == j)]), with = F],  paste0("phenotype.",j,".txt"),sep = "\t", row.names=F, quote=F, na = "NA")
    
  }

}

#--------------------------------------#
##--   create file for batch jobs   --## 
#--------------------------------------#

## create helper file to enable large-scale job arrays, careful,
## does not include male- and female-specific outcomes
tmp <- expand.grid(1:23, 1:13,  stringsAsFactors = F)
write.table(tmp, "Job.submission.REGENIE.step2.phecodes.txt", sep="\t", col.names = F, row.names = F, quote = F)

##################################################
####       prepare sex-specific outcomes      ####
##################################################

## create file for men-only outcomes
tmp <- expand.grid(1:23, 15,  stringsAsFactors = F)
write.table(tmp, "Job.submission.REGENIE.step2.men.phecodes.txt", sep="\t", col.names = F, row.names = F, quote = F)

## create file for men-only outcomes
tmp <- expand.grid(1:23, 14,  stringsAsFactors = F)
write.table(tmp, "Job.submission.REGENIE.step2.women.phecodes.txt", sep="\t", col.names = F, row.names = F, quote = F)

##################################################
####            examine job results           ####
##################################################

#--------------------------------------------#
##--       read in REGENIE log-files      --##
#--------------------------------------------#

## get all log files
log.files <- dir("../output/")
log.files <- grep("\\.log", log.files, value=T)
## all available

## loop though to get information needed
log.files <- lapply(log.files, function(x){
  
  cat("\n-----------------------------\n")
  cat(x, "\n")
  
  ## import the corresponding log file and parse information needed
  tmp         <- readLines(paste0("../output/", x))
  
  if(length(tmp) > 40){
    
    ## get the case/control numbers
    cc          <- grep("- 'date_", tmp, value=T)
    ## inflate
    cc          <- data.frame(do.call(rbind, strsplit(cc, " ")))
    ## keep only what is needed
    cc          <- cc[, c(5,6,9)]
    ## rename
    names(cc)   <- c("id", "cases", "controls")
    ## do some cleaning
    cc$id       <- gsub("'|:", "", cc$id)
    ## make numeric
    cc$cases    <- as.numeric(cc$cases)
    cc$controls <- as.numeric(cc$controls)
    
    ## search for potential warnings of fitted probabilites
    warn        <- grep("WARNING", tmp, value=T)
    # print(warn)
    ## extract the phenotypes
    warn        <- strsplit(warn, "#")
    warn        <- lapply(warn, function(k){
      ## subset to what is interesting
      k <- k[2]
      ## extract the phecode
      k <- substr(k, 1, stringr::str_locate(k, "\\)\\.")[,1]-1)
      return(k)
    })
    
    ## add indication
    cc$warn     <- cc$id %in% warn 
    
    ## add some basic characteristics
    x           <- strsplit(x, "_|\\.")[[1]]
    ## add to the data
    cc$batch    <- as.numeric(x[3])
    cc$chr      <- gsub("chr", "", x[4])
    return(cc)
    
  }else{
    cat("NO LOG-FILE for", x, "\n")
    cat("------------------------------\n")
  }
  
})
## combine into one
log.files <- do.call(rbind, log.files)
