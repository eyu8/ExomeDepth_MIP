#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(glue)
library(readr)

pd_number <- c("pd_15","pd_19","pd_22","pd_26","pd_28","pd_55","pd_65","pd_71","pd_76")
resulttable <- lapply(pd_number, function(x)readRDS(paste("raw/",x,"_result_filtered_test_v2.rds",sep="")))


args <- commandArgs(trailingOnly = TRUE)

col <- 1

for(n in 1:length(resulttable)){
        result_list <- resulttable[[n]]
    for(i in 1:length(result_list)){
        result <- result_list[[i]]
        if(class(result) == "data.frame"){
        gene <- result[result$chromosome == 6,]
        col_names <- ifelse(col == 1, TRUE, FALSE)
        col <- 0
        write_excel_csv(gene,glue("ExomeDepth_all_result_test_v2_{args[1]}.csv"),append = !col_names, col_names = col_names)
        } else {
            message(result)
        }
    }
}


