#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/CNV/ExomeDepth")


library(ExomeDepth)
library(parallel)

ncore <- 40

unrelated <- paste0("Exome.",as.character(readLines("~/runs/eyu8/QC_NeuroX/QC_cryptic_HSP_controls/HSP_unrelated.csv")),".clean.dedup.recal.bam")
hsp <- readRDS("raw/HSP_case_all_v3.rds")

ctrl.raw <- readRDS("raw/HSP_control_all_v3.rds")
ctrl.tmp <- ctrl.raw[,colnames(ctrl.raw) %in% c(colnames(ctrl.raw)[1:6],unrelated)]

ctrl_wCNV <- paste0("Exome.",as.character(readLines("HSP/HSP_controls_with_CNV_nobias")),".clean.dedup.recal.bam")

#no poor controls
ctrl.tmp2 <- ctrl.tmp
ctrl <-  ctrl.tmp2[, !(names(ctrl.tmp2) %in% ctrl_wCNV)]


### prepare the main matrix of read count data 
hsp.mat <- as.matrix(hsp[, grep(names(hsp), pattern = 'Exome*')])
ctrl.mat <- as.matrix(ctrl.tmp2[, grep(names(ctrl.tmp2), pattern = 'Exome*')])
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


nsamples <- ncol(ctrl.mat)
result <- mclapply(1:nsamples,
                   function(i){

### Create the aggregate reference set for this sample
                       ref.mat <- reference.mat[,!(colnames(reference.mat) %in% colnames(ctrl_wCNV.mat)[i])]
                       my.choice <- select.reference.set (test.counts = ctrl.mat[,i],
                                                          reference.counts = ref.mat,
                                                          data = ctrl.tmp2,
                                                          formula = 'cbind(test, reference) ~ GC',
                                                          bin.length = (ctrl.tmp2$end - ctrl.tmp2$start)/1000,
                                                          n.bins.reduced = 10000)
                       my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                                                      MAR = 1,
                                                      FUN = sum)
                       message('Now creating the ExomeDepth object')
                       all.exons <- new('ExomeDepth',
                                        test = ctrl.mat[,i],
                                        reference = my.reference.selected,
                                        data = ctrl,
                                        formula = 'cbind(test, reference) ~ GC')

############### Now call the CNVs
                       all.exons <- CallCNVs(x = all.exons,
                                             transition.probability = 10^-4,
                                             chromosome = ctrl.tmp2$space,
                                             start = ctrl.tmp2$start,
                                             end = ctrl.tmp2$end,
                                             name = ctrl.tmp2$names)
                       result <- all.exons@CNV.calls

                       if(nrow(result)> 0){
                           result$Sample = colnames(ctrl.tmp2)[i+6]
                           result$correlation = cor(all.exons@test, all.exons@reference)
                           result$numRef = length(my.choice$reference.choice)
                           result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')
                       }
                       return(result)
                   }, mc.cores=ncore)
saveRDS(result,"raw/HSP_control_result_v1_Exome.rds")






