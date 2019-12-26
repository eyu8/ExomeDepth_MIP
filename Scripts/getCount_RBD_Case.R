#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")


library(ExomeDepth)
library(tidyverse)
library(glue)

bed <- "original.targets_only.sorted.bed"
fasta <- "Reference/human_g1k_v37.fasta" 
#change loop and comment vector to use
rbd <- c("rbd_25","rbd_26","rbd_28","rbd_55","rbd_71")


for(name in rbd){

        name.temp <- str_replace(name,"_",".PD_MIP")
        case <- scan(glue('../MIP/rbd_noGBAP1_mip/{name.temp}.bam.noGBAP1.files'),character(),quote="")
        count <- getBamCounts(bed.file = bed, bam.files = case, include.chr = FALSE, referenceFasta = fasta)
        saveRDS(count,glue('raw/{name}_all.rds'))
}



