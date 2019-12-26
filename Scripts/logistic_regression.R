#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(plyr)

result <- data.frame(Gene = character(), CNV = double(), Freq_PD = double(), Freq_Ctrl = double(), pvalue = double(), correction = double(), OR = double(), OR_5 =  double(), OR_95 = double(),stringsAsFactors=FALSE)

table <- readRDS("raw/logistic_regression_table.rds")

table[table$Sex==0,]$Sex <- NA
table[table$Age==0,]$Age <- NA
table$Pheno <- table$Pheno-1
#deletion analysis

for(i in 5:length(names(table))){
    table[table[,i]==3,i] <- 0
    table[table[,i]==2,i] <- 0
}



for(i in 3:length(names(table))){
    table[,i] <- factor(table[,i])
}


for(i in 5:length(names(table))){
    counts<- ddply(table, .(table$Pheno, table[,i]), nrow)
    names(counts) <- c("Pheno","deletion")
    Freq_PD <- counts[4,3]/counts[3,3]
    Freq_Ctrl <- counts[2,3]/counts[1,3]
    gene <- colnames(table)[i]
    if(length(levels(table[,gene]))< 2){
        
        result[i-4,] = list(gene,1,NA,NA,NA,NA,NA,NA,NA)
        next
    }
    mylogit <- glm(Pheno ~ get(gene) + Age + Sex, data = table, family = binomial() )

    p <- summary(mylogit)$coefficient[2,4]
    correction <- p * 50
    or <- exp(summary(mylogit)$coefficient[2,1])
    cfd <- exp(confint(mylogit,level = 0.90))
    or5 <- cfd[2,1]
    or95 <- cfd[2,2]
    result[i-4,] = list(gene,1,Freq_PD,Freq_Ctrl,p,correction,or,or5,or95)
}

result[result$correction < 0.05,]

#duplication

for(i in 5:length(names(table))){
    table[table[,i]==1,i] <- 0
    table[table[,i]==2,i] <- 0
}



for(i in 3:length(names(table))){
    table[,i] <- factor(table[,i])
}


for(i in 5:length(names(table))){
    counts<- ddply(table, .(table$Pheno, table[,i]), nrow)
    names(counts) <- c("Pheno","duplication")
    if(nrow(counts)==4){
    Freq_PD <- counts[4,3]/counts[3,3]
    Freq_Ctrl <- counts[2,3]/counts[1,3]
    } else {
        Freq_PD <- NA
        Freq_Ctrl <- NA
    }
    gene <- colnames(table)[i]    
    if(length(levels(table[,gene]))< 2){
        result[i-4,] = list(gene,3,NA,NA,NA,NA,NA,NA,NA)
        next
    }        
    mylogit <- glm(Pheno ~ get(gene) + Age + Sex, data = table, family = binomial() )
   
    p <- summary(mylogit)$coefficient[2,4]
    correction <- p * 50
    or <- exp(summary(mylogit)$coefficient[2,1])
    cfd <- exp(confint(mylogit,level = 0.90))
    or5 <- cfd[2,1]
    or95 <- cfd[2,2]
    result[i-4,] = list(gene,3,Freq_PD,Freq_Ctrl,p,correction,or,or5,or95)
}

result[result$correction < 0.05,]

