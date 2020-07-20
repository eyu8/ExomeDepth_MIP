# ExomeDepth for Molecular Inversion Probes (MIPs) in the *PRKN* gene

## Published papers

 Papers published with release: 
 - [ExomeDepth for MIPs v1.0.0](https://github.com/gan-orlab/ExomeDepth_MIP/releases/tag/v1.0.0)
	> *Analysis of heterozygous PRKN variants and copy number variations in Parkinson’s disease.* Under review. 2020. Preprint available on medRxiv: DOI: https://doi.org/10.1101/2020.05.07.20072728
	> *Analysis of dominant and recessive parkinsonism genes in REM sleep behavior disorder.* Under review. 2020. Preprint available on medRxiv: DOI: https://doi.org/10.1101/2020.03.17.20032664


## WORKFLOW 
#### [1. Extract reads from MIPs](#1)
#### [2. Preparing MIPs filter files](#2)
#### [3. Filtering MIPs](#3)
#### [4. Call CNVs](#4)
#### [5. Annotate result](#5)

<a id="1"></a>
## 1. Extract reads from MIPs 
These 3 packages are necessary to run this CNV calling pipeline. 
```R
library(ExomeDepth)
library(tidyverse)
library(parallel)
library(readr)
```
First step is to extract reads from MIPs bam files using **ExomeDepth::getBamCounts** function. In this example, extraction was done in batches. 
```R
#bed file containing all probes
bed <- "MIPs_probes.bed"
#reference file
fasta <- "/home/eyu8/Reference/human_g1k_v37.fasta"

#list of pd batches
pd <- c("pd_15","pd_19","pd_22","pd_26","pd_28","pd_55","pd_65","pd_71","pd_76","pd_77","pd_78","pd_79")

for(name in pd){
        name.temp <- str_replace(name,"_",".PD_MIP")
        #location of PD bam files
        case <- scan(paste0('../MIP/pd_noGBAP1_mip/', name.temp, '.bam.noGBAP1.files'),character(),quote="")
        #get read count for each bam
        count <- getBamCounts(bed.file = bed, bam.files = case, include.chr = FALSE, referenceFasta = fasta)
        #save in a R binary file
        saveRDS(count,paste0('raw/', name, '_all.rds'))
}
```
<a id="2"></a>
## 2. Preparing MIPs filter files
Three filters were applied to MIPs before CNV calling:
- Exclude duplicates
- Exclude probes with low coverage
- Exclude probes from genes with low coverage
- Exclude samples with average coverage of less than 50x
```R
#list of duplicates as vector
duplicates <- readLines("../../MIP/txt/data/NeuroX_duplicate_v3.data")
duplicates <- paste0('MIP.',duplicates,'.clean.bam',sep='')

#list of probe with low coverage (<100x) as vector
badProbes <- readLines("../../MIP/txt/cov/bad_probe.cov")

#list of gene with good coverage (100x in at least 85% of the gene) as vector collapsed with "|"
#to be used as a pattern for grep in the next step
good_gene <- paste(readLines("../../MIP/txt/cov/good_gene_above_90_vps35.cov"),collapse="|")

```
<a id="3"></a>
## 3. Filtering MIPs
```R
#load read count of all PD samples
pd_number <- c("pd_15","pd_19","pd_22","pd_26","pd_28","pd_55","pd_65","pd_71","pd_76")
pd <- lapply(pd_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
#select gene that pass QC
pd_PARK2 <- lapply(pd, function(x) x[grep(rownames(x), pattern = good_gene),])
pd.list <- lapply(pd_PARK2, function(x) as(x[, colnames(x)], 'data.frame'))
#remove duplicates
pd.tmp1 <- lapply(pd.list, function(x) x[ , !(names(x) %in% duplicates)])
#remove probes that don't pass QC
pd.tmp2 <- lapply(pd.tmp1, function(x) x[!(x$name %in% badProbes),])
#remove samples with less than 50X
pd.dafr <- lapply(pd.tmp2, function(x) ifelse(length(which(colMeans(x[,-c(1:6)]) < 50)) > 0,return(x[ , -(as.vector(which(colMeans(x[,-c(1:6)]) < 50))+6)]),return(x)))


#load read count of all control samples
ctrl_number <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")
ctrl <- lapply(ctrl_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
#select gene that pass QC
ctrl_PARK2 <- lapply(ctrl, function(x) x[grep(rownames(x), pattern = good_gene),])
ctrl.list <- lapply(ctrl_PARK2, function(x) as(x[, colnames(x)], 'data.frame'))
#remove duplicates
ctrl.tmp1 <- lapply(ctrl.list, function(x) x[ , !(names(x) %in% duplicates)])
#remove probes that don't pass QC
ctrl.tmp2 <- lapply(ctrl.tmp1, function(x) x[!(x$name %in% badProbes),])
#remove samples with less than 50X
ctrl.dafr <- lapply(ctrl.tmp2, function(x) ifelse(length(which(colMeans(x[,-c(1:6)]) < 50)) > 0,return(x[ , -(as.vector(which(colMeans(x[,-c(1:6)]) < 50))+6)]),return(x)))

### prepare the main matrix of read count data
pd.mat <- as.matrix(pd.dafr[, grep(names(pd.dafr), pattern = '*.clean.bam')])
control.list.mat <- lapply(ctrl.dafr, function(x) as.matrix(x[, grep(names(x), pattern = '*.clean.bam')]))
reference.mat <- Reduce(cbind,lapply(ctrl.dafr, function(x) as.matrix(x[, grep(names(x), pattern = '*.clean.bam')])))

```
<a id="4"></a>
## 4. Call CNVs

```R
#declare number of CPU core
ncore <- 40


## start looping over each sample
for(n in 1:length(pd_number)){
    pd_count.dafr <- pd.dafr[[n]]
    pd_count.mat <- pd.mat[[n]]
    nsamples <- ncol(pd_count.mat)


    result <- mclapply(1:nsamples,
                       function(i){
### Create the aggregate reference set for this sample
                           my.choice <- select.reference.set (test.counts = pd_count.mat[,i],
                                                              reference.counts = reference.mat,
                                                              data = pd_count.dafr,
                                                              formula = 'cbind(test, reference) ~ GC')
                           my.reference.selected <- apply(X = reference.mat[, my.choice$reference.choice, drop = FALSE],
                                                          MAR = 1,
                                                          FUN = sum)
                           message('Now creating the ExomeDepth object')
                           all.exons <- new('ExomeDepth',
                                            test = pd_count.mat[,i],
                                            reference = my.reference.selected,
                                            data = pd_count.dafr,
                                            formula = 'cbind(test, reference) ~ GC')
############### Now call the CNVs
                           all.exons <- CallCNVs(x = all.exons,
                                                 transition.probability = 10^-6,
                                                 chromosome = pd_count.dafr$space,
                                                 start = pd_count.dafr$start,
                                                 end = pd_count.dafr$end,
                                                 name = pd_count.dafr$names,
                                                 expected.CNV.length = 100000)
                           result <- all.exons@CNV.calls
                           if(nrow(result)> 0){
                               result$Sample = colnames(pd_count.dafr)[i+6]
                               result$correlation = cor(all.exons@test, all.exons@reference)
                               result$numRef = length(my.choice$reference.choice)
                               result$ReferenceSet = paste(my.choice$reference.choice,collapse = ',')

                           }
                           return(result)
                       }, mc.cores=ncore)
    saveRDS(result,paste0('raw/', pd_number[n], '_result.rds'))
}
```
<a id="5"></a>
## 5. Annotate result for *PRKN*
Load the files saved in the previous step and combine them.
```R
pd_number <- c("pd_15","pd_19","pd_22","pd_26","pd_28","pd_55","pd_65","pd_71","pd_76")
resulttable <- lapply(pd_number, function(x)readRDS(paste("raw/",x,"_result.rds",sep="")))
pd.cnv <- Reduce(rbind,unlist(resulttable,recursive=FALSE))

ctrl_number <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")
resulttable <- lapply(ctrl_number, function(x)readRDS(paste("raw/",x,"_result.rds",sep="")))
control.cnv <- Reduce(rbind,unlist(resulttable,recursive=FALSE))
```
Select *PRKN* CNV
```R
#select gene of interest (start and end chosen based on probe target region)
pd.csv <- pd.cnv[pd.cnv$chromosome == 6 & pd.cnv$start > 161768590 & pd.cnv$end < 163148834,]
control.csv <- control.cnv[control.cnv$chromosome == 6 & control.cnv$start > 161768590  & control.cnv$end < 163148834,]
```
Highlight CNVs with correlation of less than 0.97
```R
#Filter column created with Bad, Cor, Keep where Bad means CNV exist in control and Cor means correlation is less than 0.97
pd.csv$Filter <- as.character("Keep")
pd.csv[pd.csv$correlation <= 0.97,]$Filter <- as.character("Cor")
```
Save the results
```R
write_excel_csv(pd.csv, "PRKN_PD.csv")
write_excel_csv(control.csv, "PRKN_control.csv")
```
