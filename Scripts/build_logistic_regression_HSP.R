#! /usr/bin/env Rscript

packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/CNV/ExomeDepth")

library(ExomeDepth)
library(readr)


args <- commandArgs(trailingOnly = TRUE)

newExons <- read.table("HSP/HSP_RefSeq_genes",as.is = T)

#poor_region <- readLines("HSP_poor_region")
#hsp.genes <- readLines("HSP_RefSeq_genes_only")[-as.numeric(poor_region)]
#hsp.genes <- readLines("HSP_RefSeq_genes_only")
resulttable <- readRDS(args[1])

hsp.cnv <- Reduce(rbind,resulttable)

resulttable <- readRDS(args[2])
#resulttable <- readRDS("raw/HSP_controlCNV.rds")
control.cnv <- Reduce(rbind,resulttable)

hsp.cnv$start.gene <- "Any"
hsp.cnv$end.gene <- "Any"
control.cnv$start.gene <- "Any"
control.cnv$end.gene <- "Any"

hsp.cnv$Filter <- as.character("Keep")
if(any(hsp.cnv$correlation <= 0.97)){hsp.cnv[hsp.cnv$correlation <= 0.97,]$Filter <- as.character("Cor")}
for(i in 1:nrow(hsp.cnv)){
    hsp.cnv$start.gene[i] <- newExons[as.numeric(hsp.cnv$start.p[i]),]$V4
    hsp.cnv$end.gene[i] <- newExons[as.numeric(hsp.cnv$end.p[i]),]$V4
} 

for(i in 1:nrow(control.cnv)){
    control.cnv$start.gene[i] <- newExons[as.numeric(control.cnv$start.p[i]),]$V4
    control.cnv$end.gene[i] <- newExons[as.numeric(control.cnv$end.p[i]),]$V4
}

for(i in 1:nrow(control.cnv)){
    remove <- hsp.cnv[,1] == control.cnv[i,1] & hsp.cnv[,2] == control.cnv[i,2] &  hsp.cnv[,3] == control.cnv[i,3]
    if(any(remove)){hsp.cnv[remove,]$Filter <- as.character("Bad")}
}

hsp.cnv$S_Number <- substr(hsp.cnv$Sample,7,12)
control.cnv$S_Number <- substr(control.cnv$Sample,7,12)

write_excel_csv(hsp.cnv,args[3])
write_excel_csv(control.cnv,args[4])

