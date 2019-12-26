#! /usr/bin/env Rscript


packrat::init("~/runs/eyu8/library/ExomeDepth")
setwd("~/runs/eyu8/data/ExomeDepth")

library(glue)
library(readr)

result_list <- readRDS("raw/control_result_filtered_test_v2.rds")


args <- commandArgs(trailingOnly = TRUE)

col <- 1
for(i in 1:length(result_list)){
    result <- result_list[[i]]
    if(class(result) == "data.frame"){
        gene <- result[result$chromosome == 6,]
        col_names <- ifelse(col == 1, TRUE, FALSE)
        col <- 0
        write_excel_csv(gene,glue("ExomeDepth_control_test_v2_{args[1]}.csv"),append = !col_names, col_names = col_names)
    } else {
        message(result)
    }
}
