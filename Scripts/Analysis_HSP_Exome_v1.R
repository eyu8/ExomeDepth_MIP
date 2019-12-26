#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")


library(ExomeDepth)
library(parallel)

ncore <- 40

unrelated <- paste0("Exome.",as.character(readLines("~/runs/eyu8/QC_cryptic_HSP_controls/HSP_unrelated.csv")),".clean.dedup.recal.bam")
hsp <- readRDS("raw/HSP_case_all_v3.rds")

ctrl.raw <- readRDS("raw/HSP_control_all_v3.rds")
ctrl.tmp <- ctrl.raw[,colnames(ctrl.raw) %in% c(colnames(ctrl.raw)[1:6],unrelated)]
ctrl <- ctrl.tmp

### prepare the main matrix of read count data 
hsp.mat <- as.matrix(hsp[, grep(names(hsp), pattern = 'Exome*')])
reference.mat <- as.matrix(ctrl[, grep(names(ctrl), pattern = 'Exome*')])


## start looping over each sample

nsamples <- ncol(hsp.mat)
result <- mclapply(1:nsamples,
                   function(i){

### Create the aggregate reference set for this sample
                       my.choice <- select.reference.set (test.counts = hsp.mat[,i],
                                                          reference.counts = reference.mat,
                                                          data = hsp,
                                                          formula = 'cbind(test, reference) ~ GC',
                                                          bin.length = (hsp$end - hsp$start)/1000,
                                                          n.bins.reduced = 10000)
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
                                             transition.probability = 10^-4,
                                             chromosome = hsp$space,
                                             start = hsp$start,
                                             end = hsp$end,
                                             name = hsp$names)
                       result <- all.exons@CNV.calls

                       if(nrow(result)> 0){
                           result$Sample = colnames(hsp)[i+6]
                           result$correlation = cor(all.exons@test, all.exons@reference)
                           result$numRef = length(my.choice$reference.choice)
                           result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')
                       }
                       return(result)
                   }, mc.cores=ncore)
saveRDS(result,"raw/HSP_result_v1_Exome.rds")


nsamples <- ncol(reference.mat)
result <- mclapply(1:nsamples,
                   function(i){

### Create the aggregate reference set for this sample
                       my.choice <- select.reference.set (test.counts = reference.mat[,i],
                                                          reference.counts = reference.mat[,-i],
                                                          data = ctrl,
                                                          formula = 'cbind(test, reference) ~ GC',
                                                          bin.length = (ctrl$end - ctrl$start)/1000,
                                                          n.bins.reduced = 10000)
                       my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                                                      MAR = 1,
                                                      FUN = sum)
                       message('Now creating the ExomeDepth object')
                       all.exons <- new('ExomeDepth',
                                        test = reference.mat[,i],
                                        reference = my.reference.selected,
                                        data = ctrl,
                                        formula = 'cbind(test, reference) ~ GC')

############### Now call the CNVs
                       all.exons <- CallCNVs(x = all.exons,
                                             transition.probability = 10^-4,
                                             chromosome = ctrl$space,
                                             start = ctrl$start,
                                             end = ctrl$end,
                                             name = ctrl$names)
                       result <- all.exons@CNV.calls

                       if(nrow(result)> 0){
                           result$Sample = colnames(ctrl)[i+6]
                           result$correlation = cor(all.exons@test, all.exons@reference)
                           result$numRef = length(my.choice$reference.choice)
                           result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')
                       }
                       return(result)
                   }, mc.cores=ncore)
saveRDS(result,"raw/HSP_control_result_v1_Exome.rds")






