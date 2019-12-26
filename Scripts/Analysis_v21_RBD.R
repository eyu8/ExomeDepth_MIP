#! /usr/bin/env Rscript

packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(ExomeDepth)
library(tidyverse)
library(glue)
library(parallel)

ncore <- 40

duplicates <- readLines("../MIP/txt/data/NeuroX_duplicate_v3.data")
duplicates <- paste0('MIP.',duplicates,'.clean.bam',sep='')

badProbes <- readLines("../MIP/txt/cov/bad_probe.cov")
#poorSamples <- readLines("../MIP/txt/files/PARK2_bad.csv")
#poorSamples <- c(poorSamples,"S33024")
#poorSamples <- paste0('MIP.',poorSamples,'.clean.bam',sep='')

good_gene <- paste(readLines("../MIP/txt/cov/good_gene_above_90_vps35.cov"),collapse="|")

#sample_less_50 <- readLines("../MIP/txt/cov/all_gene_less_50.cov")

pd_number <- c("rbd_25","rbd_26","rbd_28","rbd_55","rbd_71")
pd <- lapply(pd_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
pd_PARK2 <- lapply(pd, function(x) x[grep(rownames(x), pattern = good_gene),])
pd.list <- lapply(pd_PARK2, function(x) as(x[, colnames(x)], 'data.frame'))
pd.tmp1 <- lapply(pd.list, function(x) x[ , !(names(x) %in% duplicates)])
pd.tmp2 <- lapply(pd.tmp1, function(x) x[!(x$name %in% badProbes),])
pd.dafr <- lapply(pd.tmp2, function(x) ifelse(length(which(colMeans(x[,-c(1:6)]) < 50)) > 0,return(x[ , -(as.vector(which(colMeans(x[,-c(1:6)]) < 50))+6)]),return(x)))
#pd.tmp2 <- lapply(pd.tmp1, function(x) x[ , !(names(x) %in% sample_less_50)])
#pd.dafr <- lapply(pd.tmp1, function(x) x[ , !(names(x) %in% poorSamples)])


ctrl_number <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")
ctrl <- lapply(ctrl_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
ctrl_PARK2 <- lapply(ctrl, function(x) x[grep(rownames(x), pattern = good_gene),])
ctrl.list <- lapply(ctrl_PARK2, function(x) as(x[, colnames(x)], 'data.frame'))
ctrl.tmp1 <- lapply(ctrl.list, function(x) x[ , !(names(x) %in% duplicates)])
ctrl.tmp2 <- lapply(ctrl.tmp1, function(x) x[!(x$name %in% badProbes),])
ctrl.dafr <- lapply(ctrl.tmp2, function(x) ifelse(length(which(colMeans(x[,-c(1:6)]) < 50)) > 0,return(x[ , -(as.vector(which(colMeans(x[,-c(1:6)]) < 50))+6)]),return(x)))

#ctrl.tmp2 <- lapply(ctrl.tmp1, function(x) x[ , !(names(x) %in% sample_less_50)])
#ctrl.dafr <- lapply(ctrl.tmp1, function(x) x[ , !(names(x) %in% poorSamples)])

### prepare the main matrix of read count data
pd.mat <- lapply(pd.dafr, function(x) as.matrix(x[, grep(names(x), pattern = '*.clean.bam')]))
control.list.mat <- lapply(ctrl.dafr, function(x) as.matrix(x[, grep(names(x), pattern = '*.clean.bam')]))
reference.mat <- Reduce(cbind,lapply(ctrl.dafr, function(x) as.matrix(x[, grep(names(x), pattern = '*.clean.bam')])))
#reference.mat <- reference.mat[,-length(colnames(reference.mat))]



## start looping over each sample
for(n in 1:length(pd_number)){
    pd_count.dafr <- pd.dafr[[n]]
    pd_count.mat <- pd.mat[[n]]
    nsamples <- ncol(pd_count.mat)

   
    result <- mclapply(1:nsamples,
                       function(i){
### Create the aggregate reference set for this sample
                           my.choice <- select.reference.set (test.counts = pd_count.mat[,i],
                                                              reference.counts = reference.mat,
                                                              data = pd_count.dafr,
                                                              formula = 'cbind(test, reference) ~ GC')
                           my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                                                          MAR = 1,
                                                          FUN = sum)
                           message('Now creating the ExomeDepth object')
                           all.exons <- new('ExomeDepth',
                                            test = pd_count.mat[,i],
                                            reference = my.reference.selected,
                                            data = pd_count.dafr,
                                            formula = 'cbind(test, reference) ~ GC')
############### Now call the CNVs
                           all.exons <- CallCNVs(x = all.exons,
                                                 transition.probability = 10^-6,
                                                 chromosome = pd_count.dafr$space,
                                                 start = pd_count.dafr$start,
                                                 end = pd_count.dafr$end,
                                                 name = pd_count.dafr$names,
                                                 expected.CNV.length = 100000)
                           result <- all.exons@CNV.calls
                           if(nrow(result)> 0){
                               result$Sample = colnames(pd_count.dafr)[i+6]
                               result$correlation = cor(all.exons@test, all.exons@reference)
                               result$numRef = length(my.choice$reference.choice)
                               result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')

                           }
                           return(result)
                       }, mc.cores=ncore)
    saveRDS(result,glue('raw/{pd_number[n]}_result_test_v21.rds'))
}


for(n in 1:length(ctrl_number)){
    candidate_control.dafr <- ctrl.dafr[[n]]
    control.mat <- control.list.mat[[n]]
    nsamples <- ncol(control.mat)
    result <- mclapply(1:nsamples,
                       function(i){
#### Create the aggregate reference set for this sample
                           looRef.mat <- reference.mat[,colnames(reference.mat) != colnames(control.mat)[i]]
                           my.choice <- select.reference.set (test.counts = control.mat[,i],
                                                              reference.counts = looRef.mat,
                                                              data = candidate_control.dafr,
                                                              formula = 'cbind(test, reference) ~ GC')
                           my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                                                          MAR = 1,
                                                          FUN = sum)
                           message('Now creating the ExomeDepth object')
                           all.exons <- new('ExomeDepth',
                                            test = control.mat[,i],
                                            reference = my.reference.selected,
                                            data = candidate_control.dafr,
                                            formula = 'cbind(test, reference) ~ GC')
################ Now call the CNVs
                           all.exons <- CallCNVs(x = all.exons,
                                                 transition.probability = 10^-6,
                                                 chromosome = candidate_control.dafr$space,
                                                 start = candidate_control.dafr$start,
                                                 end = candidate_control.dafr$end,
                                                 name = candidate_control.dafr$names,
                                                 expected.CNV.length = 100000)
                           result <- all.exons@CNV.calls
                           if(nrow(result) > 0){
                               result$Sample = colnames(candidate_control.dafr)[i+6]
                               result$correlation = cor(all.exons@test, all.exons@reference)
                               result$numRef = length(my.choice$reference.choice)
                               result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')
                           }
                           return(result)
                       }, mc.cores=ncore)
    saveRDS(result,glue('raw/{ctrl_number[n]}_result_test_v21.rds'))
}
