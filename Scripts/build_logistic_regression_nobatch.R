#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/CNV/ExomeDepth")

library(ExomeDepth)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
#load results from ExomeDepth
resulttable <- readRDS(args[1])
pd.cnv <- Reduce(rbind,resulttable)

#select gene of interest (start and end chosen based on probe target region)
pd.csv <- pd.cnv[pd.cnv$chromosome == as.numeric(args[2]) & pd.cnv$start > as.numeric(args[3]) & pd.cnv$end < as.numeric(args[4]),]
#Filter column created with Bad, Cor, Keep where Bad means CNV exist in control and Cor means correlation is less than 0.97
pd.csv$Filter <- as.character("Keep")
if(any(pd.csv$correlation <= 0.97)){pd.csv[pd.csv$correlation <= 0.97,]$Filter <- as.character("Cor")}
pd.csv$S_Number <- substr(pd.csv$Sample,5,10)

write_excel_csv(pd.csv, args[5])

