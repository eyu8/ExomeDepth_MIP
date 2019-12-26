#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(ExomeDepth)
#library(glue)
library(parallel)

ncore <- 40


ROI <- as.numeric(readLines("HSP_RefSeq_row"))
unrelated <- paste0("Exome.",as.character(readLines("~/runs/eyu8/QC_cryptic_HSP_controls/HSP_unrelated.csv")),".clean.dedup.recal.bam")
hsp.raw <- readRDS("raw/HSP_case_all.rds")
hsp.tmp <- hsp.raw[ROI,]

ctrl.raw <- readRDS("raw/HSP_control_all.rds")
ctrl.tmp <- ctrl.raw[ROI,colnames(ctrl.raw) %in% c(colnames(ctrl.raw)[1:6],unrelated)]

#poor_region <- unique(as.vector(which(rowMeans(hsp.tmp[,-c(1:6)]) < 15)), as.vector(which(rowMeans(ctrl.tmp[,-c(1:6)]) < 15)))
#no poor ctrls
#poor_ctrls <- as.vector(which(colMeans(ctrl.tmp[,-c(1:6)]) < 15))+6
#ctrl_wCNV <- readLines("HSP_controls_with_CNV")


hsp <- hsp.tmp
ctrl <- ctrl.tmp
#ctrl <- ctrl.tmp[-poor_region,-poor_ctrls]
#ctrl <-  ctrl.tmp2[, !(names(ctrl.tmp2) %in% ctrl_wCNV)]

### prepare the main matrix of read count data
hsp.mat <- as.matrix(hsp[, grep(names(hsp), pattern = 'Exome*')])
reference.mat <- as.matrix(ctrl[, grep(names(ctrl), pattern = 'Exome*')])
#ctrl_wCNV.mat <- as.matrix(ctrl.tmp2[, grep(names(ctrl.tmp2), pattern = 'Exome*')])

## start looping over each sample

nsamples <- ncol(hsp.mat)
result <- mclapply(1:nsamples,
                   function(i){

### Create the aggregate reference set for this sample
                       my.choice <- select.reference.set (test.counts = hsp.mat[,i],
                                                          reference.counts = reference.mat,                      
                                                          data = hsp,
                                                          formula = 'cbind(test, reference) ~ GC')
                       my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                                                      MAR = 1,
                                                      FUN = sum)
                       message('Now creating the ExomeDepth object')
                       all.exons <- new('ExomeDepth',
                                        test = hsp.mat[,i],
                                        reference = my.reference.selected,
                                        data = hsp,
                                        formula = 'cbind(test, reference) ~ GC')
############### Now call the CNVs
                       all.exons <- CallCNVs(x = all.exons,
                                             transition.probability = 10^-2,
                                             chromosome = hsp$space,
                                             start = hsp$start,
                                             end = hsp$end,
                                             name = hsp$names,
                                             expected.CNV.length = 100000)
                       result <- all.exons@CNV.calls

                       if(nrow(result)> 0){
                           result$Sample = colnames(hsp)[i+6]
                           result$correlation = cor(all.exons@test, all.exons@reference)
                           result$numRef = length(my.choice$reference.choice)
                           result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')                           
                       }
                       return(result)
                   }, mc.cores=ncore)
saveRDS(result,"raw/HSP_result_v4_all.rds")


nsamples <- ncol(reference.mat)
result <- mclapply(1:nsamples,
                   function(i){
                           #### Create the aggregate reference set for this sample
                       my.choice <- select.reference.set (test.counts = reference.mat[,i],
                                                          reference.counts = reference.mat[,-i],
                                                          data = ctrl,
                                                          formula = 'cbind(test, reference) ~ GC')
                       my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                                                      MAR = 1,
                                                      FUN = sum)
                       message('Now creating the ExomeDepth object')
                       all.exons <- new('ExomeDepth',
                                            test = reference.mat[,i],
                                            reference = my.reference.selected,
                                            data = ctrl,
                                            formula = 'cbind(test, reference) ~ GC')
                       ################ Now call the CNVs
                       all.exons <- CallCNVs(x = all.exons,
                                             transition.probability = 10^-2,
                                             chromosome = ctrl$space,
                                             start = ctrl$start,
                                             end = ctrl$end,
                                             name = ctrl$names,
                                             expected.CNV.length = 100000)
                       result <- all.exons@CNV.calls
                       if(nrow(result) > 0){
                           result$Sample = colnames(ctrl)[i+6]
                           result$correlation = cor(all.exons@test, all.exons@reference)
                           result$numRef = length(my.choice$reference.choice)
                           result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')
                       }
                       return(result)
                   }, mc.cores=ncore)
saveRDS(result,'raw/HSP_control_result_v4_all.rds')



