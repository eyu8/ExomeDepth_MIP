#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")


args <- commandArgs(trailingOnly = TRUE)

start <- args[1]
end <- args[2]


library(ExomeDepth)
data(exons.hg19)

#bed <- "RefSeq.Genes.v37.all_exons.no_flanks.uniq_regions.sort.bed"
fasta <- "Reference/human_g1k_v37.fasta" 
#change loop and comment vector to use


case <- scan("../hsp/hsp.bam",character(),quote="")
count <- getBamCounts(bed.frame = exons.hg19, bam.files =  case[start:end], include.chr = FALSE, referenceFasta = fasta)
saveRDS(count,paste0('raw/HSP_case/HSP_case_v3_',start,'.rds'))
