#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")


library(ExomeDepth)
library(tidyverse)
library(glue)

#bed file containing all probes
bed <- "original.targets_only.sorted.bed"
#reference file
fasta <- "/home/eyu8/Reference/human_g1k_v37.fasta" 
#list of pd batches
#pd <- c("pd_15","pd_19","pd_22","pd_26","pd_28","pd_55","pd_65","pd_71","pd_76")
pd <- c("pd_78")

for(name in pd){

        name.temp <- str_replace(name,"_",".PD_MIP")
        #location of PD bam files
        case <- scan(glue('../MIP/pd_noGBAP1_mip/{name.temp}.bam.noGBAP1.files'),character(),quote="")
        #get read count for each bam
        count <- getBamCounts(bed.file = bed, bam.files = case, include.chr = FALSE, referenceFasta = fasta)
        #save in a R binary file
        saveRDS(count,glue('raw/{name}_all.rds'))
}



