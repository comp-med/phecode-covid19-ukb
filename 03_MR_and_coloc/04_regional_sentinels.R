####################################################
#### collate regional lead signals for COVID-19 ####
#### Maik Pietzner                   04/09/2023 ####
####################################################

rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(igraph)
require(doMC)
require(igraph)

#########################################
####    import regional sentinels    ####
#########################################

#------------------------------#
##--    COVID-19 outcomes   --##
#------------------------------#

## get all regional sentinel files
ii           <- dir("../regional_results/")
## restrict to files needed
ii           <- grep("[0-9]_regional", ii, value=T)

## get the file header
hd           <- fread(cmd="zcat ../../01_stats/A2_ALL_exUKBB.tsv.gz | head -1")
## rename one
names(hd)[1] <- "CHR"

## import the results
res.regional <- lapply(ii, function(x){
  
  ## import
  tmp        <- fread(paste0("../regional_results/", x))
  print(head(tmp))
  ## add header
  names(tmp) <- c(names(hd), "region_start", "region_end")
  ## add endpoint
  tmp[, endpoint := gsub("_regional_sentinels.txt", "", x)]
  ## return
  return(tmp)
  
})
## combine everything
res.regional <- do.call(rbind, res.regional)

#------------------------------#
##--        Long COVID      --##
#------------------------------#

## import results for LongCOVID
tmp.regional        <- fread("../regional_results/st_cas_vs_bd_con_regional_sentinels.txt")
## adopt names to match to COVID-19 stats
names(tmp.regional) <- c("rsid", "CHR", "POS", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "REF", "ALT", "all_meta_AF", "all_inv_var_meta_p", "all_inv_var_meta_effective", "region_start", "region_end")
## endpoint
tmp.regional[, endpoint := "LC"]

## lift over position to HG19 (manual via rsid and )
tmp.regional$rsid
tmp.regional[, POS := c(43153864, 41483390)]

## combine everything
res.regional        <- plyr::rbind.fill(res.regional, tmp.regional)
res.regional        <- as.data.table(res.regional)
## create MarkerName
res.regional[, MarkerName := paste0("chr", CHR, ":", POS, "_", pmin(REF, ALT), "_", pmax(REF,ALT))]

#------------------------------#
##--    map to UKB stats    --##
#------------------------------#  

## import UKB SNPs to merge
ukb.snps     <- fread("<path to file>")

## merge
res.regional <- merge(res.regional, ukb.snps[, c("MarkerName", "rsid")], by="MarkerName", all.x=T, suffix=c(".hgi", ".ukb"))
## delete what is no longer needed
rm(ukb.snps); gc(reset=T)
## one SNP missing
res.regional[ is.na(rsid.ukb) ]

## export for query in phecode GWAS stats
write.table(unique(na.omit(res.regional$rsid.ukb)), "phecode.variants.lookup.covid.20230912.txt", sep="\t", col.names = F, row.names = F, quote = F)
write.table(unique(na.omit(res.regional$rsid.hgi)), "covid.variants.lookup.covid.20230921.txt", sep="\t", col.names = F, row.names = F, quote = F)

#########################################
####          perform coloc          ####
#########################################

#-----------------------------#
##--    process look-up    --##
#-----------------------------#

## import look-up files (n=81 out of 82 results)
ii                  <- dir("../tmpdir/")
ii                  <- grep("allchr.results.gz", ii, value=T)
## do in parallele
registerDoMC(20)
## import
covid.lookup        <- mclapply(ii, function(x){
  ## import
  tmp    <- fread(paste0("../tmpdir/", x))
  ## add phenotype
  tmp$id <- gsub("phecode\\.|\\.allchr.results.gz.variant.lookup.results", "", x)
  ## return results
  return(tmp)
}, mc.cores = 20)
## combine
covid.lookup        <- do.call(rbind, covid.lookup)
## get effect estimates
names(covid.lookup) <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA", "id")
## create MarkerName
covid.lookup[, MarkerName := paste0("chr", as.numeric(gsub("X", 23, CHROM)), ":", as.numeric(GENPOS), "_", pmin(ALLELE0,ALLELE1), "_", pmax(ALLELE0,ALLELE1))]

## how many findings with at least suggestive evidence
nrow(covid.lookup[ LOG10P > 6])
## n = 173

## add covid-19 statistics (A2, B2 and long COVID, only)
covid.lookup        <- merge(covid.lookup, res.regional[ endpoint %in% c("A2", "B2", "LC"), c("MarkerName", "endpoint", "REF", "ALT", "all_meta_AF", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_inv_var_meta_p")],
                             by="MarkerName") 
length(unique(covid.lookup$MarkerName)) ## n=64
length(unique(res.regional[ endpoint %in% c("A2", "B2", "LC"), ]$MarkerName)) ## n=65

## align effect direction
covid.lookup[, all_inv_var_meta_beta := ifelse(ALLELE1 == ALT, all_inv_var_meta_beta, -all_inv_var_meta_beta)]
## adjust allele frequency accordingly
covid.lookup[, all_meta_AF := ifelse(ALLELE1 == ALT, all_meta_AF, 1-all_meta_AF)]
## drop what is no longer needed
covid.lookup$ALT <- covid.lookup$REF <- NULL

## overall alignment of effect estimates
plot(BETA ~ all_inv_var_meta_beta, covid.lookup[LOG10P > 5])
## consider heatmap for outcomes that provide evidence for coloc

#-----------------------------#
##--      input coloc      --##
#-----------------------------#

## write to file (at least suggestive evidence)
## omit MHC region! - 6:25500000-34000000
write.table(covid.lookup[LOG10P > 5 & !(CHROM == 6 & GENPOS >= 25500000 & GENPOS <= 34000000), c("id", "endpoint", "CHROM", "GENPOS")], "input.phecode.coloc.txt", sep="\t", row.names = F, col.names = F, quote = F)

#-----------------------------#
##--     import results    --##
#-----------------------------#

## import results
ii <- dir("../output/")
## reduce to output of interest
ii <- grep("susie", ii, value=T)
## seems to be all there

## import results
res.coloc <- lapply(ii, function(x) fread(paste0("../output/", x)))
## combine
res.coloc <- do.call(plyr::rbind.fill, res.coloc)
res.coloc <- as.data.table(res.coloc)
## filter 
res.coloc <- res.coloc[ ld.sens^2 > .8 & ld.top > .8 & PP.H4.abf > .8 ]

#-----------------------------#
##-- merge COVID-19 clumps --##
#-----------------------------#

## file of SNPs to be queried
write.table(unique(res.regional$rsid.ukb), "../tmpdir/candidate.snps.txt", sep="\t", row.names = F, col.names = F, quote = F)

## run scripts to extract from autosomes
system(paste("./../scripts/obtain_snps.sh", "covid"))

## import SNPs
snps         <- fread("../tmpdir/snp.covid.dosage.transpose")
## import SNP information
snps.info    <- fread("../tmpdir/snp.covid.info") ## n=78 found (X-chromosome missing, ignore for now)
## create markername 
snps.info[, MarkerName := paste0("chr", chromosome, "_", position, "_", pmin(alleleA, alleleB), "_", pmax(alleleA, alleleB))]
## assign names to SNP file
names(snps)  <- c("IID", snps.info$MarkerName)
## import inclusion list
in.list      <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/genotypes/sample_inclusion/qctool_pass_EUR_panukbb_all_unrelated.incl")
## subset
snps         <- snps[ IID %in% in.list$ID_1]

## compute LD-matrix
ld.mat       <- cor(snps[, -1])^2
## convert to data table
ld.mat[ lower.tri(ld.mat, diag = F)] <- NA
ld.mat       <- as.data.frame(as.table(ld.mat))
## create network
ld.sub       <- graph_from_data_frame(subset(ld.mat, Freq >= .7))
## get all separate components
ld.sub       <- components(ld.sub)$membership
## convert to data frame
ld.sub       <- data.table(ID=names(ld.sub), R2.group=ld.sub)
## adopt names
ld.sub[, MarkerName := sub("_", ":", ID)]

## add to the data
res.regional <- merge(res.regional, ld.sub[, c("MarkerName", "R2.group")], by="MarkerName", all.x=T)
table(res.regional$R2.group, useNA="always")

## delete what is no longer needed
rm(snps); rm(snps.info); gc(reset=T)

#-----------------------------#
##--     add locus data    --##
#-----------------------------#

## add locus data to coloc results
res.coloc    <- merge(res.coloc, unique(res.regional[, c("MarkerName", "R2.group")]), by.x="covid19.lead", by.y="MarkerName")
table(res.coloc$R2.group, useNA = "always")
#   4    8   16   24   27   42   47   56   60 <NA> 
#   4    2    1    5   20    2    6    2   12    0 

## get most likely gene assignments for each locus 
## group from the flagship paper
write.table(unique(res.coloc[, c("covid19.lead", "R2.group")]), "R2.group.COVID19.assignment.20230921.txt", sep="\t", row.names=F)
## import annotation
r2.anno      <- fread("R2.group.COVID19.assignment.20230921.txt")
## add to coloc results
res.coloc    <- merge(res.coloc, unique(r2.anno[, c("R2.group", "Gene.assignment")]), by="R2.group") 

## create effects based on the COVID-19 risk increasing allele
res.coloc[, BETA.phecode.increasing := sign(BETA.covid)*BETA.phecode]
res.coloc[, BETA.covid.increasing := abs(BETA.covid)]

#-----------------------------#
##--    export network     --##
#-----------------------------#

## --> create edges for a bipartite graph <-- ##
tmp1        <- res.coloc[, c("id.phecode", "Gene.assignment", "BETA.phecode.increasing"), with=F]
tmp2        <- res.coloc[, c("id.covid", "Gene.assignment", "BETA.covid.increasing"), with=F]
## edit names
names(tmp1) <- names(tmp2) <- c("id", "Gene.assignment", "BETA")
e.tmp       <- rbind(tmp1, tmp2)
## add another column to ease plotting
e.tmp$lty   <- ifelse(e.tmp$BETA > 0, 1, 2)
## make unique
e.tmp       <- unique(e.tmp)

## --> add vertex label <-- ##
tmp1        <- unique(res.coloc[, c("id.covid", "phenotype.covid")])
## add category
tmp1[, category := "infectious diseases"]
## same for phecode
tmp2        <- unique(res.coloc[, c("id.phecode", "phenotype.phecode")])
tmp2        <- merge(tmp2, unique(res.phecode[, c("id", "category")]), by.x = "id.phecode", by.y = "id")
## temp. combine
names(tmp1) <- names(tmp2) <- c("id", "label", "category")
tmp         <- rbind(tmp1, tmp2)
tmp$type    <- "phenotype"
## add colour
tmp         <- merge(tmp, unique(res.phecode[, c("category", "cl")]))
## gene label
tmp1        <- unique(res.coloc[, c("Gene.assignment", "R2.group")])
tmp1[, rsid := sapply(R2.group, function(x) res.regional$rsid.hgi[which(res.regional$R2.group == x)][1])]
## keep only what is needed and rename
tmp1[, label := paste0(rsid, " (", Gene.assignment, ")")]
tmp1[, type := "locus"]
tmp1[, cl := "white"]
tmp1[, category := "locus"]
tmp1        <- tmp1[, c("category", "Gene.assignment", "label", "type", "cl")]
## edit names
names(tmp1) <- c("category", "id", "label", "type", "cl")
## combine
v.tmp       <- rbind(tmp, tmp1)

## --> write to file <-- ##
write.table(v.tmp, "Nodes.COVID19.phecode.network.20230921.txt", sep="\t", row.names = F, quote = F)
write.table(e.tmp, "Edges.COVID19.phecode.network.20230921.txt", sep="\t", row.names = F, quote = F)

#-----------------------------#
##--  create forest plot   --##
#-----------------------------#

## --> phecode stats <-- ##

## prepare data to do so
tmp1        <- unique(res.coloc[, c("Gene.assignment", "R2.group")])
tmp1[, rsid := sapply(R2.group, function(x){
  ## get all
  ii <- res.regional$rsid.ukb[which(res.regional$R2.group == x)]
  ## select the one included in the UKB data
  return(ii[ ii %in% covid.lookup$ID][1])
})]
## phecodes
tmp1        <- merge(tmp1, unique(res.coloc[, c("id.phecode", "Gene.assignment")]), by="Gene.assignment")
tmp1        <- merge(tmp1, unique(covid.lookup[, c("id", "MarkerName", "ID", "ALLELE0", "ALLELE1", "BETA", "SE", "LOG10P")]), by.x=c("id.phecode", "rsid"), by.y=c("id", "ID"), all.x=T)

## --> COVID-19 stats <-- ##

## import look-up
tmp2 <- fread("../output/covid19.covid19.lookup.results")
names(tmp2)   <- c("misc", names(res.regional)[4:22])
## create chromosome and phenotype
tmp2$CHROM    <- as.numeric(gsub(".*:", "", tmp2$misc)) 
tmp2$endpoint <- gsub("<path to file>/|_ALL_exUKBB.tsv.gz", "", gsub(":.*", "", tmp2$misc))
## create MarkerName
tmp2[, MarkerName := paste0("chr", as.numeric(gsub("X", 23, CHROM)), ":", as.numeric(POS), "_", pmin(REF, ALT), "_", pmax(REF, ALT))]
## add effect alleles for phecode data
tmp2          <- merge(tmp2, unique(covid.lookup[, c("MarkerName", "ID", "ALLELE0", "ALLELE1")]), by="MarkerName")
## add locus 
tmp2          <- merge(tmp2, tmp1[, c("MarkerName", "R2.group", "Gene.assignment")], by="MarkerName")
## add columns to merge with network nodes
tmp2[, BETA := ifelse( ALT == ALLELE1, all_inv_var_meta_beta, -all_inv_var_meta_beta)]
tmp2[, LOG10P := -log10(all_inv_var_meta_p)]
## keep only what is needed
tmp2          <- tmp2[ endpoint %in% c("A2", "B2"), c("endpoint", "ID", "Gene.assignment", "R2.group", "MarkerName", "ALLELE0", "ALLELE1", "BETA", "all_inv_var_meta_sebeta", "LOG10P")]
## rename
names(tmp2)   <- c("id.phecode", "rsid", "Gene.assignment", "R2.group", "MarkerName", "ALLELE0", "ALLELE1", "BETA", "SE", "LOG10P")
## subset to rsids included
tmp2          <- unique(tmp2)

## --> combine <-- ##
res.forest <- plyr::rbind.fill(tmp1, tmp2)
res.forest <- as.data.table(res.forest)
## order for plotting
res.forest <- res.forest[ order(Gene.assignment, id.phecode),]
## drop ABO und FUT2 for now
res.forest <- res.forest[ !(Gene.assignment %in% c("ABO", "FUT2"))]
res.forest[, plt.srt := 1:nrow(res.forest)]

## create a column for COVID-19 increasing effects
tmp.bet    <- unique(res.forest[ id.phecode == "A2", c("rsid", "BETA")])
tmp.bet[, increase := sign(BETA)]
## merge again
res.forest <- merge(res.forest, tmp.bet[, c("rsid", "increase")])
## BETA for the plot
res.forest[, BETA.plot := BETA * increase]
## reorder
res.forest <- res.forest[order(Gene.assignment, id.phecode), ]
  
# png("../graphics/Forest.plot.locus.effects.COVID19.20230921.png", width = 8, height = 8, res = 600, units = "cm")
pdf("../graphics/Forest.plot.locus.effects.COVID19.20230921.pdf", width = 3.15, height = 3.15)
## plotting parameters
par(mar=c(1.5,7,.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, mgp=c(.6,0,0), xaxs="i", yaxs="i")

## create tmp data frame to ease plotting
tmp <- do.call(data.frame, aggregate(plt.srt ~ rsid, res.forest, function(x) c(min(x), mean(x), max(x))))
tmp <- tmp[ order(tmp$plt.srt.1),]

## empty plot
plot(c(-1.2, 1), c(.5, nrow(res.forest)+.5), type = "n", ylab="", yaxt="n", xaxt="n",
     xlab="Odds ratio (95%-CI)", ylim=rev(c(.5, nrow(res.forest)+.5)))
## add axis
axis(1, at=log(c(.5,.75,1,1.5,2,3)), labels = c(.5,.75,1,1.5,2,3), lwd = .5)
## add rectangle to divide
pm <- par("usr")
rect(pm[1], tmp$plt.srt.1-.5, pm[2], tmp$plt.srt.3+.5, col=c("white", "grey80"), border = NA)
## abline no effect
abline(v=0, lwd=.5)

## add confidence intervals
arrows(res.forest$BETA.plot - 1.96 * res.forest$SE, 1:nrow(res.forest), res.forest$BETA.plot + 1.96 * res.forest$SE, 1:nrow(res.forest),
       length = 0, lwd = .5, col=ifelse(res.forest$id.phecode %in% c("A2", "B2"), "#00A4CC", "#F95700"))
## add point estimates
points(res.forest$BETA.plot, 1:nrow(res.forest), pch=22, lwd=.3, bg=ifelse(res.forest$id.phecode %in% c("A2", "B2"), "#00A4CC", "#F95700"),
       col="grey20", cex=.7)

## add label
text(pm[1], 1:nrow(res.forest), cex=.35, xpd=NA, pos=2, offset = .1, 
     labels = sapply(res.forest$id.phecode, function(x) v.tmp$label[which(v.tmp$id == x)]))

## add gene label
text(pm[1], tmp$plt.srt.1+.5, pos=4, labels = sapply(tmp$rsid, function(x) paste0(x, "\n(",res.forest$Gene.assignment[which(res.forest$rsid == x)][1], ")")),
     cex=.35)

## legend
legend("bottomright", bty="n", cex=.5, lty=0, pch=22, pt.lwd=.5, pt.bg=c("#00A4CC", "#F95700"),
       legend = c("COVID-19", "Disease UKB"), pt.cex=.8)

## close device
dev.off()

#-----------------------------#
##-- create heatmap across --##
#-----------------------------#

## add labels to ease some downstream analysis
res.forest <- merge(res.forest, v.tmp, by.x="id.phecode", by.y="id")

## --> get covid stats for all COVID19 hits <-- ##

## import look-up
tmp2          <- fread("../output/covid19.covid19.lookup.results")
names(tmp2)   <- c("misc", names(res.regional)[4:22])
## create chromosome and phenotype
tmp2$CHROM    <- as.numeric(gsub(".*:", "", tmp2$misc)) 
tmp2$endpoint <- gsub("<path to file>/|_ALL_exUKBB.tsv.gz", "", gsub(":.*", "", tmp2$misc))
## create MarkerName
tmp2[, MarkerName := paste0("chr", as.numeric(gsub("X", 23, CHROM)), ":", as.numeric(POS), "_", pmin(REF, ALT), "_", pmax(REF, ALT))]
## add effect alleles for phecode data
tmp2          <- merge(tmp2, unique(covid.lookup[, c("MarkerName", "ID", "ALLELE0", "ALLELE1")]), by="MarkerName")
## add locus 
tmp2          <- merge(tmp2, unique(res.regional[, c("MarkerName", "R2.group")]), by="MarkerName")
## align effect estimates to the effect allele in UKB
tmp2[, BETA := ifelse( ALT == ALLELE1, all_inv_var_meta_beta, -all_inv_var_meta_beta)]
tmp2[, LOG10P := -log10(all_inv_var_meta_p)]
## keep only what is needed
tmp2          <- tmp2[ endpoint %in% c("A2", "B2"), c("endpoint", "ID", "R2.group", "MarkerName", "ALLELE0", "ALLELE1", "BETA", "all_inv_var_meta_sebeta", "LOG10P")]
## rename
names(tmp2)   <- c("id", "ID", "R2.group", "MarkerName", "ALLELE0", "ALLELE1", "BETA", "SE", "LOG10P")
## creaete indicator what to drop
tmp2[, ind := 1:.N, by=c("R2.group")]
## keep only what is needed (6 outcomes per variant)
tmp2          <- tmp2[ ind <= 2 ]
## n=49 variants, one with only one outcome

## delete what is no longer needed
tmp2[, ind := NULL]

## --> get phecode stats for all COVID19 hits <-- ##

## limit to selected endpoints (A2, B2, date_480.2, date_495, date_496, date_502, date_504, date_714.1)
res.heat   <- covid.lookup[ id %in% c("date_480.2", "date_495", "date_496", "date_502", "date_504", "date_714.1"), c("id", "MarkerName", "ID", "ALLELE0", "ALLELE1", "BETA", "SE", "LOG10P")]
res.heat   <- unique(res.heat)
## SNPs also in the COVID-19 data
res.heat   <- res.heat[ ID %in% tmp2$ID]
## add R2 groups to filter
res.heat   <- merge(res.heat, unique(res.regional[, c("MarkerName", "R2.group")]), by = "MarkerName")
## n = 49 variants

## --> combine into one data frame <-- ##

## combine 
res.heat   <- rbind(tmp2, res.heat)

## create a column for COVID-19 increasing effects
tmp.bet    <- unique(res.heat[ id == "A2", c("ID", "BETA")])
tmp.bet[, increase := sign(BETA)]
## merge again
res.heat   <- merge(res.heat, tmp.bet[, c("ID", "increase")])
## BETA for the plot
res.heat[, BETA.plot := BETA * increase]
## delete what is no longer needed
res.heat[, increase := NULL]
## reshape to enable clustering
res.heat   <- reshape(res.heat, idvar = c("MarkerName", "ID", "ALLELE0", "ALLELE1", "R2.group"), timevar = "id", direction = "wide")

## --> plot simple heatmap <-- ##

## add chromosome and position to ease plotting
res.heat            <- merge(unique(covid.lookup[, c("MarkerName", "CHROM", "GENPOS")]), res.heat)
## order accordingly
res.heat            <- res.heat[ order(CHROM, GENPOS),]

## start device
pdf("../graphics/Heatmap.locus.effects.COVID19.20230921.pdf", width = 6.3, height = 3.15)
# png("../graphics/Heatmap.locus.effects.COVID19.20230921.png", width = 16, height = 8, res = 600, units = "cm")
## plotting parameters
par(mar=c(5,8,.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, mgp=c(.6,0,0), xaxs="i", yaxs="i", lwd=.5)

## define colour gradient
col.grad <- colorRampPalette(c("#00A4CC", "white", "#F95700"))(230)

## empty plot
plot(c(.5, nrow(res.heat)+.5), c(.5, length(var.clust)+2.5), type="n", xlab="", ylab = "", yaxt = "n", xaxt =  "n")

## add coloured rectangles
k <- 0
for(j in c("A2", "B2", var.clust)){
  ## increase counter
  k <- k+1
  ## add coloured rectangle
  rect(1:nrow(res.heat)-.5, k-.5, 1:nrow(res.heat)+.5, k+.5, border="white", lwd=.2, col = col.grad[110 + ceiling(unlist(res.heat[, paste0("BETA.plot.", j), with=F])*100)])
}

## surround if gw-significant
k <- 0
for(j in c("A2", "B2", var.clust)){
  ## increase counter
  k <- k+1
  ## go through each entry
  for(jj in 1:nrow(res.heat)){
    if(!is.na(res.heat[jj, paste0("LOG10P.", j), with=F]) & res.heat[jj, paste0("LOG10P.", j), with=F] > 7.3){
      ## add coloured rectangle
      rect(jj-.5, k-.5, jj+.5, k+.5, border="grey20", lwd=.2, col = col.grad[110 + ceiling(unlist(res.heat[jj, paste0("BETA.plot.", j), with=F])*100)])
    }
  }
}

## get plotting coordinates
pm <- par("usr")

## add traits
text(pm[1], 1:length(c("A2", "B2", var.clust)), cex=.5, pos=2, xpd=NA, 
     labels = sapply(c("A2", "B2", var.clust), function(x) v.tmp$label[which(v.tmp$id == x)]))

## add genetic variants
text(1:nrow(res.heat), pm[3]-(pm[4]-pm[3])*.02, cex=.5, pos=2, xpd=NA, srt=90, labels = res.heat$ID, offset=0)

## --> add color legend <-- ##

## length
ll <- seq(pm[1]-(pm[2]-pm[1])*.3, pm[1]-(pm[2]-pm[1])*.1, length.out = length(col.grad))
## rectangle for the colours
rect(ll-(ll[2]-ll[1])/2, pm[3]-(pm[4]-pm[3])*.2, ll+(ll[2]-ll[1])/2, pm[3]-(pm[4]-pm[3])*.1, border=NA, col=col.grad, xpd=NA)
## box
rect(ll[1]-(ll[2]-ll[1])/2, pm[3]-(pm[4]-pm[3])*.2, ll[length(ll)]+(ll[2]-ll[1])/2, pm[3]-(pm[4]-pm[3])*.1, border="black", col=NA, lwd=.3, xpd=NA)
## add header
text(pm[1]-(pm[2]-pm[1])*.3, pm[3]-(pm[4]-pm[3])*.08, cex=.6, labels = "Effect estimate", pos=4,
     offset = .2, xpd=NA)
## simple axis
text(ll[c(10, 60, 120, 180, 220)], pm[3]-(pm[4]-pm[3])*.2,
     labels=c(-1,-.5, 0, .5, 1), pos=1, cex=.4, offset = .1, xpd=NA)

## close device
dev.off()

#-----------------------------#
##--      Suppl. Table     --##
#-----------------------------#

## write results to file
write.table(res.coloc, "Results.COVID19.loci.phecode.coloc.20231002.txt", sep="\t", row.names=F)
