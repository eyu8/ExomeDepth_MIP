#!/bin/bash

module load gcc/7.3.0 r/3.5.2
Rscript Scripts/build_logistic_regression_nobatch.R $1 6 161768456 163148994 result_v21/PARK2_$2.csv 

Rscript Scripts/build_logistic_regression_nobatch.R $1 4 90645125 90759594 result_v21/SNCA_$2.csv 

