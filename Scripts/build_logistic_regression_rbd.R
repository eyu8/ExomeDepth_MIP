#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(ExomeDepth)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

pd_number <- c("rbd_25","rbd_26","rbd_28","rbd_55","rbd_71")
#pd_number <- c("pd_55")
resulttable <- lapply(pd_number, function(x)readRDS(paste("raw/",x,"_result_test_v",args[1],".rds",sep="")))

pd.cnv <- Reduce(rbind,unlist(resulttable,recursive=FALSE))

#ctrl_number <- c("control_55")
ctrl_number <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")
resulttable <- lapply(ctrl_number, function(x)readRDS(paste("raw/",x,"_result_test_v",args[1],".rds",sep="")))
control.cnv <- Reduce(rbind,unlist(resulttable,recursive=FALSE))
pd.csv <- pd.cnv[pd.cnv$chromosome == as.numeric(args[2]) & pd.cnv$start > as.numeric(args[3]) & pd.cnv$end < as.numeric(args[4]),]
pd.csv$Filter <- as.character("Keep")
control.csv <- control.cnv[control.cnv$chromosome == as.numeric(args[2]) & control.cnv$start > as.numeric(args[3]) & control.cnv$end < as.numeric(args[4]),]

if(any(pd.csv$correlation <= 0.97)){
    pd.csv[pd.csv$correlation <= 0.97,]$Filter <- as.character("Cor")
}

for(i in 1:nrow(control.csv)){
    remove <- pd.csv[,1] == control.csv[i,1] & pd.csv[,2] == control.csv[i,2] &  pd.csv[,3] == control.csv[i,3]
    if(any(remove)){pd.csv[remove,]$Filter <- as.character("Bad")}
}

pd.csv$S_Number <- substr(pd.csv$Sample,5,10)
control.csv$S_Number <- substr(control.csv$Sample,5,10)

write_excel_csv(pd.csv, args[5])
write_excel_csv(control.csv, args[6])

