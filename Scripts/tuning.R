#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/CNV/ExomeDepth")

library(ExomeDepth)
library(readr)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

expected.CNV.length <- args[1]
#10^-6
transition.probability <- args[2]
version <- args[3]
core <- args[4]

duplicates <- readLines("../../MIP/txt/data/NeuroX_duplicate_v3.data")
duplicates <- paste0('MIP.',duplicates,'.clean.bam',sep='')
badProbes <- readLines("../../MIP/txt/cov/bad_probe.cov")
#badProbes <- c()

#poorSamples <- readLines("../MIP/txt/files/PARK2_bad.csv")
#poorSamples <- paste0('MIP.',poorSamples,'.clean.bam',sep='')

MLPA_cnv <- paste0('MIP.',readLines("txt/MLPA_cnv"),'.clean.bam',sep='')
MLPA_wild <- paste0('MIP.',readLines("txt/MLPA_wild"),'.clean.bam',sep='')

good_gene <- paste(readLines("../../MIP/txt/cov/good_gene_above_90_vps35.cov"),collapse="|")

ctrl_number <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")
ctrl <- lapply(ctrl_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
ctrl_PARK2 <- lapply(ctrl, function(x)x[grep(rownames(x), pattern = good_gene),])
ctrl.list <- lapply(ctrl_PARK2, function(x) as(x[, colnames(x)], 'data.frame'))
ctrl.merged.bam <- Reduce(cbind,lapply(ctrl.list, function(x) as(x[, colnames(x)], 'data.frame')))
ctrl.merged.dafr <- ctrl.merged.bam[, c(1:6,grep(names(ctrl.merged.bam), pattern = '*.clean.bam'))]
ctrl.tmp1 <- ctrl.merged.dafr[ , !(names(ctrl.merged.dafr) %in% duplicates)]
ctrl.tmp2 <- ctrl.tmp1[!(ctrl.tmp1$names %in% badProbes),]
ctrl.dafr <- ctrl.tmp2[ , -(as.vector(which(colMeans(ctrl.tmp2[,-c(1:6)]) < 50))+6)]
ctrl.bam <- ctrl.dafr[, grep(names(ctrl.dafr), pattern = '*.clean.bam')]
ctrl.ref <- ctrl.dafr[ , !(names(ctrl.dafr) %in% c(MLPA_cnv,MLPA_wild))]



pd_number <- c("pd_15","pd_19","pd_22","pd_26","pd_28","pd_55","pd_65","pd_71","pd_76")
pd <- lapply(pd_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
pd_PARK2 <- lapply(pd, function(x) x[grep(rownames(x), pattern = good_gene),])
pd.tmp <- Reduce(cbind,lapply(pd_PARK2, function(x) as(x[, colnames(x)], 'data.frame')))
pd.list <-  pd.tmp[!(pd.tmp[,5] %in% badProbes),]
pd.bam <- pd.list[, c(1:6,grep(names(pd.list), pattern = '*.clean.bam'))]
all.bam <- cbind(pd.bam,ctrl.bam)
CNV_samples.dafr <- all.bam[ , (names(all.bam) %in% c(MLPA_cnv,colnames(all.bam)[1:6]))]
noCNV_samples.dafr <- all.bam[ , (names(all.bam) %in% c(MLPA_wild,colnames(all.bam)[1:6]))]
SNCA.dafr <- all.bam[ , (names(all.bam) %in% c("MIP.S33485.clean.bam","MIP.S33460.clean.bam",colnames(all.bam)[1:6]))]
### prepare the main matrix of read count data
CNV_samples.mat <- as.matrix(CNV_samples.dafr[, grep(names(CNV_samples.dafr), pattern = '*.clean.bam')])
noCNV_samples.mat <- as.matrix(noCNV_samples.dafr[, grep(names(noCNV_samples.dafr), pattern = '*.clean.bam')])
reference.mat <- as.matrix(ctrl.ref[, grep(names(ctrl.ref), pattern = '*.clean.bam')])
SNCA.mat <- noCNV_samples.mat[,colnames(noCNV_samples.mat) %in% c("MIP.S33485.clean.bam","MIP.S33460.clean.bam")]

### start looping over each sample
getResult <- function(mat, dafr, nsamples,exp,tp,core){
    mclapply(1:nsamples, function(i){
#### Create the aggregate reference set for this sample
    message("Sample number: ",colnames(mat)[i])
    my.choice <- select.reference.set (test.counts = mat[,i],
                                       reference.counts = reference.mat,
                                       data = dafr,
                                       formula = 'cbind(test, reference) ~ GC')

    my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                                   MAR = 1,
                                   FUN = sum)
    
    message('Now creating the ExomeDepth object')
    all.exons <- new('ExomeDepth',
                     test = mat[,i],
                     reference = my.reference.selected,
                     data = dafr,
                     formula = 'cbind(test, reference) ~ GC')
################ Now call the CNVs

    message("Reference for this samples:",length(my.choice$reference.choice))
   
    all.exons <- CallCNVs(x = all.exons,
                          transition.probability = as.numeric(tp),
                          chromosome = dafr$space,
                          start = dafr$start,
                          end = dafr$end,
                          name = dafr$names,
                          expected.CNV.length = as.numeric(exp))
    result <- all.exons@CNV.calls
    if(nrow(result)> 0){
        result$Sample = colnames(dafr)[i+6]
        result$correlation = cor(all.exons@test, all.exons@reference)
        result$numRef = length(my.choice$reference.choice)
        result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')        
    }
   if(colnames(mat)[i] == "MIP.S33485.clean.bam"){
   pdf(paste0(colnames(mat)[i],".pdf"))
   plot(all.exons,
          #sequence = '6',
          sequence = '4',
          #sequence = '1',
          xlim = c(90645250 , 90759447),
          #xlim = c(161768590 , 163148834 ),
          #xlim = c(20959948, 20978004),
          count.threshold = 20,
          main = 'SNCA deletion',
          #main = 'PINK1',
          #main = 'PARK2 gene',
          cex.lab = 0.8,
          with.gene = TRUE)
    dev.off()
   }
    return(result)
}, mc.cores=core)}


#CNV_result_MLPA <- as.data.frame(read_csv("txt/CNV_result_MLPA.csv"))
#resulttable <- getResult(CNV_samples.mat,CNV_samples.dafr,ncol(CNV_samples.mat),expected.CNV.length,transition.probability,core)
################################ 
#score_list <- lapply(resulttable, function(result){
#                         if(nrow(result) == 0){
#                             c(0,1)
#                         } else {
#                             MLPA <- CNV_result_MLPA[CNV_result_MLPA$Sample == substr(result$Sample,5,10)[1],]
#                             gene <- result[result$chromosome==MLPA$chr & result$start > MLPA$start & result$end < MLPA$end,]
#                             TP<- ifelse(nrow(gene) > 0, 1, 0)
#                             FN <- 1 - TP
#                             c(TP,FN)
#                         }
#})

#s <- Reduce('+',score_list)

test <- getResult(SNCA.mat,SNCA.dafr,ncol(SNCA.mat),expected.CNV.length,transition.probability,core)
MLPA_genes <- as.data.frame(read_csv("txt/PARK2_genes.csv"))
#resulttable <- getResult(noCNV_samples.mat,noCNV_samples.dafr,ncol(noCNV_samples.mat),expected.CNV.length,transition.probability,core)
score_list <- lapply(resulttable, function(result){
                         if(nrow(result) == 0){
                             c(1,0)
                         } else {
                             test <- 0
                             for(i in 1:nrow(result)){
                                 cnv <- result[i,]
                                 checkFP <- MLPA_genes$chr == cnv$chromosome & MLPA_genes$start < cnv$start & MLPA_genes$end > cnv$end
                                 if(checkFP){
                                     test <- test + 1
                                 }
                             }
                             FP <- test
                             TN <- 1 - FP
                             c(TN,FP)
                         }
})

score <- Reduce('+',score_list)
F1 <- 2*s[1]/(2*s[1] + s[2] + score[2])
FPR <- score[2]/(score[1]+score[2])
message("F1=",F1)
message("Accuracy=",s[1]/(s[1]+s[2]))
message("FPR=",FPR)
#message("Score=",score[2]," ",score[1])
