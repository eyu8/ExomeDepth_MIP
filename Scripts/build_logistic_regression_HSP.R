#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(ExomeDepth)
library(readr)
data(exons.hg19)


args <- commandArgs(trailingOnly = TRUE)

newExons <- exons.hg19[order(as.numeric(exons.hg19$chromosome),exons.hg19$start,exons.hg19$end),]
row.names(newExons) <- 1:nrow(newExons)

#poor_region <- readLines("HSP_poor_region")
#hsp.genes <- readLines("HSP_RefSeq_genes_only")[-as.numeric(poor_region)]
#hsp.genes <- readLines("HSP_RefSeq_genes_only")
resulttable <- readRDS(args[1])

hsp.cnv <- Reduce(rbind,resulttable)

resulttable <- readRDS(args[2])
control.cnv <- Reduce(rbind,resulttable)

hsp.cnv$start.gene <- "Any"
hsp.cnv$end.gene <- "Any"
control.cnv$start.gene <- "Any"
control.cnv$end.gene <- "Any"

hsp.cnv$Filter <- as.character("Keep")
if(any(hsp.cnv$correlation <= 0.97)){hsp.cnv[hsp.cnv$correlation <= 0.97,]$Filter <- as.character("Cor")}
for(i in 1:nrow(hsp.cnv)){
    hsp.cnv$start.gene[i] <- newExons[as.numeric(hsp.cnv$start.p[i]),]$name
    hsp.cnv$end.gene[i] <- newExons[as.numeric(hsp.cnv$end.p[i]),]$name
} 

for(i in 1:nrow(control.cnv)){
    control.cnv$start.gene[i] <- newExons[as.numeric(control.cnv$start.p[i]),]$name
    control.cnv$end.gene[i] <- newExons[as.numeric(control.cnv$end.p[i]),]$name
}

for(i in 1:nrow(control.cnv)){
    remove <- hsp.cnv[,1] == control.cnv[i,1] & hsp.cnv[,2] == control.cnv[i,2] &  hsp.cnv[,3] == control.cnv[i,3]
    if(any(remove)){hsp.cnv[remove,]$Filter <- as.character("Bad")}
}


write_excel_csv(hsp.cnv,args[3])
write_excel_csv(control.cnv,args[4])

