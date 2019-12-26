#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(ExomeDepth)
library(tidyverse)
library(glue)
library(parallel)


pd <- "pd_71"
pd_count <- readRDS(paste("raw/",pd,".rds",sep=""))
reference <- readRDS(paste("raw/all.control.rds",sep=""))
pd_count.dafr <- as(pd_count[, colnames(pd_count)], 'data.frame')
reference.dafr <- as(reference[, colnames(reference)], 'data.frame')



pd_count.mat <- as.matrix(pd_count.dafr[, grep(names(pd_count.dafr), pattern = '*.clean.bam')])
reference.mat <- as.matrix(reference.dafr[, grep(names(reference.dafr), pattern = '*.clean.bam')])

nsamples <- ncol(pd_count.mat)

i = 493
my.choice <- select.reference.set (test.counts = pd_count.mat[,i],
                                   reference.counts = reference.mat
                                   data =  pd_count.mat[,i])

my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                               MAR = 1,
                               FUN = sum)

message('Now creating the ExomeDepth object')
all.exons <- new('ExomeDepth',
                 test = pd_count.mat[,i],
                 reference = my.reference.selected,
                 formula = 'cbind(test, reference) ~ 1')

################ Now call the CNVs

all.exons <- CallCNVs(x = all.exons,
                          transition.probability = 10^-4,
                          chromosome = pd_count.dafr$space,
                          start = pd_count.dafr$start,
                          end = pd_count.dafr$end,
                          name = pd_count.dafr$names,
                          expected.CNV.length = 100000)

all.exons@CNV.calls


