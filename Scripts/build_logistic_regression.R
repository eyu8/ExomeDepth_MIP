#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(ExomeDepth)

pd_number <- c("pd_15","pd_19","pd_22","pd_26","pd_28","pd_55","pd_65","pd_71","pd_76")
#pd_number <- c("pd_55")
resulttable <- lapply(pd_number, function(x)readRDS(paste("raw/",x,"_result_test_v7.rds",sep="")))

pd.cnv <- Reduce(rbind,unlist(resulttable,recursive=FALSE))

#ctrl_number <- c("control_55")
ctrl_number <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")
resulttable <- lapply(ctrl_number, function(x)readRDS(paste("raw/",x,"_result_test_v7.rds",sep="")))
control.cnv <- Reduce(rbind,unlist(resulttable,recursive=FALSE))



S_Number <- c( unique(pd.cnv$Sample), unique(control.cnv$Sample))

fileName <- "../MIP/txt/data/PD_case_control_age_sex_pheno.data"

MIP.age.sex.pheno <- read.delim(fileName, header=TRUE, sep=" ", stringsAsFactors=FALSE)

MIP.age.sex.pheno$S_Number <- paste0("MIP.",MIP.age.sex.pheno$S_Number,".clean.bam")

result <- MIP.age.sex.pheno[MIP.age.sex.pheno$S_Number %in% S_Number,] 


result$PARK7 <- 0
result$ATP13A2 <- 0
result$PINK1 <- 0
result$GBA <- 0
result$SYT11 <- 0
result$RAB25 <- 0
result$RAB7L1 <- 0
result$SLC41A1 <- 0
result$PM20D1 <- 0
result$SIPA1L2 <- 0
result$HTRA2 <- 0
result$TMEM163 <- 0
result$ACMSD <- 0
result$STK39 <- 0
result$GIGYF2 <- 0
result$DNAJC13 <- 0
result$MCCC1 <- 0
result$LAMP3 <- 0
result$EIF4G1 <- 0
result$GAK <- 0
result$TMEM175 <- 0
result$DGKQ <- 0
result$BST1 <- 0
result$UCHL1 <- 0
result$SCARB2 <- 0
result$FAM47E <- 0
result$SNCA <- 0
result$PARK2 <- 0
result$GPNMB <- 0
result$FGF20 <- 0
result$ITGA8 <- 0
result$PSAP <- 0
result$INPP5F <- 0
result$SMPD1 <- 0
result$LRRK2 <- 0
result$CCDC62 <- 0
result$HIP1R <- 0
result$GCH1 <- 0
result$VPS13C <- 0
result$SETD1A <- 0
result$STX1B <- 0
result$VPS35 <- 0
result$SREBF1 <- 0
result$MAPT <- 0
result$NPC1 <- 0
result$RIT2 <- 0
result$DDRGK1 <- 0
result$USP25 <- 0
result$FBXO7 <- 0
result$PLA2G6 <- 0



fileName <- "original.targets_only.sorted.bed"
bed <- read.csv(fileName,header=FALSE, sep="\t",stringsAsFactors=FALSE)
names(bed) <- c("Chr","Start","End","Gene")
for(i in 1:nrow(pd.cnv)){
    message("PD: ",i)
    target <- pd.cnv[i,]
    only_chr <- bed[bed$Chr == target$chromosome,]
    gene <- only_chr[only_chr$Start == (target$start)-1,"Gene"]
    cnv <- result[result$S_Number %in% target$Sample,gene]
    if(length(cnv) == 0){
        next
    }
    if(cnv == 0){
        if(target$type == "deletion"){
            result[result$S_Number %in% target$Sample,gene] <- 1
        } else{
            result[result$S_Number %in% target$Sample,gene] <- 3
        }
    } else if(cnv == 1){
        if(target$type == "duplication"){
            result[result$S_Number %in% target$Sample,gene] <- 2
        }
    } else if(cnv == 3){
        if(target$type == "deletion"){
            result[result$S_Number %in% target$Sample,gene] <- 2
        }
    }
}


for(i in 1:nrow(control.cnv)){
    message("control: ",i)
    target <- control.cnv[i,]
    only_chr <- bed[bed$Chr == target$chromosome,]
    gene <- only_chr[only_chr$Start == (target$start)-1,"Gene"]
    cnv <- result[result$S_Number %in% target$Sample,gene]
    if(length(cnv) == 0){
        next
    }
    if(cnv == 0){
        if(target$type == "deletion"){
            result[result$S_Number %in% target$Sample,gene] <- 1
        } else{
            result[result$S_Number %in% target$Sample,gene] <- 3
        }
    } else if(cnv == 1){
        if(target$type == "duplication"){
            result[result$S_Number %in% target$Sample,gene] <- 2
        }
    } else if(cnv == 3){
        if(target$type == "deletion"){
            result[result$S_Number %in% target$Sample,gene] <- 2
        }
    }
}

saveRDS(result,'raw/logistic_regression_table.rds')
