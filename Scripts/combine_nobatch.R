#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/CNV/ExomeDepth/Scripts")

library(ExomeDepth)

args <- commandArgs(trailingOnly = TRUE) 
n <- args[1]
prefix <- args[2]
total <- args[3]
sample_raw <- lapply(1:n, function(x)readRDS(paste(prefix,x,".rds",sep="")))
sample_list <- lapply(sample_raw, function(x) as(x[, colnames(x)], 'data.frame'))
sample <- Reduce(cbind,sample_list)
sample_seq <- sample[, grep(names(sample), pattern = '*.clean.bam')]
sample_all <- cbind(sample[,1:6], sample_seq)
rownames(sample_all) <- rownames(sample_raw[[1]])
saveRDS(sample_all,total)

