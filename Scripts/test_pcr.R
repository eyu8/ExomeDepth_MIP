#! /usr/bin/env Rscript

packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(readr)

args <- commandArgs(trailingOnly = TRUE)

MyData <- read.csv(file=args[1], header=TRUE, sep=",", stringsAsFactors=FALSE)

myGene <- read.csv(file=args[2], header=FALSE, sep="\t", stringsAsFactors=FALSE)


MyData$PCR <- 0

for(i in 1:nrow(MyData)){
    start <- MyData[i,"start"]
    end <- MyData[i,"end"]
    if(any(myGene$V2 == (start - 1)) & any(myGene$V3 == end)){
    row1 <- which(myGene$V2 == (start - 1))
    row2 <- which(myGene$V3 == end)
    MyData[i,"PCR"] <- as.numeric(myGene[row2+1,2]) - as.numeric(myGene[row1-1,3])
    }
}

write_excel_csv(MyData, args[1])
