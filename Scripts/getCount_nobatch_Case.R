#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/CNV/ExomeDepth/Scripts")


args <- commandArgs(trailingOnly = TRUE)

start <- args[2]
end <- args[3]


library(ExomeDepth)

bed <- "../original.targets_only.sorted.bed"
fasta <- "../Reference/human_g1k_v37.fasta" 
#change loop and comment vector to use


case <- scan(args[1],character(),quote="")
count <- getBamCounts(bed.file = bed, bam.files =  case[start:end], include.chr = FALSE, referenceFasta = fasta)
saveRDS(count,paste0(args[4],start,'.rds'))
