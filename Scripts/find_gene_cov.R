#! /usr/bin/env Rscript

packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/CNV/ExomeDepth")

library(ExomeDepth)

#load reads from MIPs
pd_number <- c("pd_15","pd_19","pd_22","pd_26","pd_28","pd_55","pd_65","pd_71","pd_76")
pd <- lapply(pd_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
pd_all <-  Reduce(cbind,lapply(pd,function(x) as(x[, colnames(x)], 'data.frame')))


ctrl_number <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")
ctrl <- lapply(ctrl_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
ctrl_all <-  Reduce(cbind,lapply(ctrl,function(x) as(x[, colnames(x)], 'data.frame')))

#merge case and controls
all_mips <- cbind(pd_all,ctrl_all)

#remove meta-data
mip <- all_mips[,grep(names(all_mips), pattern = "MIP*")]
#add meta-data back
final_mip <- cbind(all_mips[,1:6],mip)


#average coverage per probe
probeMeanAll <- rowMeans(final_mip[,-c(1:6)])
#standard deviation per probe
probeSDAll <- apply(final_mip[,-c(1:6)],1,sd)

#load gene of interest
gene <- read.csv("../../MIP/MIP_PD_library.bed", sep = "\t", stringsAsFactors = FALSE, header= FALSE )
names(gene) <- c("chr","start","end","gene")

#dp of interest
dp <- 100

#% coverage per gene
cov <- sapply(gene$gene,function(g){probeMean <- probeMeanAll[grep(final_mip$name, pattern = paste0(g,".*"))]; percentCov <- 100 * length(which( probeMean > dp))/length(probeMean)})

#cov <- sapply(gene$gene,function(g){probeMean <- probeMeanAll[grep(final_mip$name, pattern = paste0(g,".*"))]})


probeMean.dafr <- data.frame("probe" = final_mip$name, "mean" = as.integer(probeMeanAll), "ratio" = probeSDAll/probeMeanAll)

#probeMean.dafr[probeMean.dafr$mean < 100,]

#names(which(cov > 85))

probeMean.85.dafr <- probeMean.dafr[grep(probeMean.dafr$probe, pattern = paste(names(which(cov > 85)),collapse="|")),]

bad <- probeMean.85.dafr[probeMean.85.dafr$mean < 100,]

bad[grep(bad$probe, pattern = "PARK2|GBA|LRRK2|PINK1|SNCA|VPS35"),]
