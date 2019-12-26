#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(ExomeDepth)

hsp_raw <- lapply(1:366, function(x)readRDS(paste("raw/HSP_case/HSP_case_v3_",x,".rds",sep="")))
hsp_list <- lapply(hsp_raw, function(x) as(x[, colnames(x)], 'data.frame'))
hsp <- Reduce(cbind,hsp_list)
hsp_exome <- hsp[, grep(names(hsp), pattern = 'Exome.*')]
hsp_all <- cbind(hsp[,1:6], hsp_exome)
saveRDS(hsp_all,'raw/HSP_case_all_v3.rds')

ctrl_raw <- lapply(1:1274, function(x)readRDS(paste("raw/HSP_control/HSP_control_v3_",x,".rds",sep="")))
ctrl_list <- lapply(ctrl_raw, function(x) as(x[, colnames(x)], 'data.frame'))
ctrl <- Reduce(cbind,ctrl_list)
ctrl_exome <- ctrl[, grep(names(ctrl), pattern = 'Exome.*')]
ctrl_all <- cbind(ctrl[,1:6], ctrl_exome)

saveRDS(ctrl_all,'raw/HSP_control_all_v3.rds')
