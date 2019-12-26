#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")


library(ExomeDepth)
library(tidyverse)
library(glue)

#bed file for all probes
bed <- "original.targets_only.sorted.bed"
fasta <- "/home/eyu8/Reference/human_g1k_v37.fasta" 
ctrl <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")

for(name in ctrl){

        name.temp <- str_replace(name,"_",".PD_MIP")
        #file contianing all bam file location
        case <- scan(glue('../MIP/control_noGBAP1_mip/{name.temp}.bam.noGBAP1.files'),character(),quote="")
        count <- getBamCounts(bed.file = bed, bam.files = case, include.chr = FALSE, referenceFasta = fasta)
        saveRDS(count,glue('raw/{name}_all.rds'))
}


