#! /usr/bin/env Rscript

#load required library
packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/CNV/ExomeDepth")

library(ExomeDepth)
library(tidyverse)
library(glue)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

pd_name <- args[1]
disease = args[2]
#declare number of CPU core
ncore <- 40

#list of duplicates
duplicates <- readLines("../../MIP/txt/data/NeuroX_duplicate_v3.data")
duplicates <- paste0('MIP.',duplicates,'.clean.bam',sep='')

#list of probe with low coverage (<100x)
badProbes <- readLines("../../MIP/txt/cov/bad_probe.cov")

#list of gene with good coverage (100x in at least 85% of the gene)
good_gene <- paste(readLines("../../MIP/txt/cov/good_gene_above_90_vps35.cov"),collapse="|")

#load read count of all PD samples
pd <- readRDS(pd_name)
#select gene that pass QC
pd_PARK2 <- pd[grep(rownames(pd), pattern = good_gene),]
pd.list <- as(pd_PARK2[, colnames(pd_PARK2)], 'data.frame')
#remove duplicates
pd.tmp1 <- pd.list[ , !(names(pd.list) %in% duplicates)]
#remove probes that don't pass QC
pd.tmp2 <- pd.tmp1[!(pd.tmp1$name %in% badProbes),]
#remove samples with less than 50X
if(length(which(colMeans(pd.tmp2[,-c(1:6)]) < 50)) > 0){
    pd.dafr <- pd.tmp2[ , -(as.vector(which(colMeans(pd.tmp2[,-c(1:6)]) < 50))+6)]
}else{
    pd.dafr <- pd.tmp2
}

#load read count of all control samples
ctrl_number <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")
ctrl <- lapply(ctrl_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
#select gene that pass QC
ctrl_PARK2 <- lapply(ctrl, function(x) x[grep(rownames(x), pattern = good_gene),])
ctrl.list <- lapply(ctrl_PARK2, function(x) as(x[, colnames(x)], 'data.frame'))
#remove duplicates
ctrl.tmp1 <- lapply(ctrl.list, function(x) x[ , !(names(x) %in% duplicates)])
#remove probes that don't pass QC
ctrl.tmp2 <- lapply(ctrl.tmp1, function(x) x[!(x$name %in% badProbes),])
#remove samples with less than 50X
ctrl.dafr <- lapply(ctrl.tmp2, function(x) ifelse(length(which(colMeans(x[,-c(1:6)]) < 50)) > 0,return(x[ , -(as.vector(which(colMeans(x[,-c(1:6)]) < 50))+6)]),return(x)))


### prepare the main matrix of read count data
pd.mat <- as.matrix(pd.dafr[, grep(names(pd.dafr), pattern = '*.clean.bam')])
control.list.mat <- lapply(ctrl.dafr, function(x) as.matrix(x[, grep(names(x), pattern = '*.clean.bam')]))
reference.mat <- Reduce(cbind,lapply(ctrl.dafr, function(x) as.matrix(x[, grep(names(x), pattern = '*.clean.bam')])))



## start looping over each sample
    nsamples <- ncol(pd.mat)
   
    result <- mclapply(1:nsamples,
                       function(i){
### Create the aggregate reference set for this sample
                           my.choice <- select.reference.set (test.counts = pd.mat[,i],
                                                              reference.counts = reference.mat,
                                                              data = pd.dafr,
                                                              formula = 'cbind(test, reference) ~ GC')
                           my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                                                          MAR = 1,
                                                          FUN = sum)
                           message('Now creating the ExomeDepth object')
                           all.exons <- new('ExomeDepth',
                                            test = pd.mat[,i],
                                            reference = my.reference.selected,
                                            data = pd.dafr,
                                            formula = 'cbind(test, reference) ~ GC')
############### Now call the CNVs
                           all.exons <- CallCNVs(x = all.exons,
                                                 transition.probability = 10^-6,
                                                 chromosome = pd.dafr$space,
                                                 start = pd.dafr$start,
                                                 end = pd.dafr$end,
                                                 name = pd.dafr$names,
                                                 expected.CNV.length = 100000)
                           result <- all.exons@CNV.calls
                           if(nrow(result)> 0){
                               result$Sample = colnames(pd.dafr)[i+6]
                               result$correlation = cor(all.exons@test, all.exons@reference)
                               result$numRef = length(my.choice$reference.choice)
                               result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')

                           }
                           return(result)
                       }, mc.cores=ncore)
    saveRDS(result,glue('raw/{disease}_result.rds'))
