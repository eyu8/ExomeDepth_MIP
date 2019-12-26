#!/bin/bash

Rscript build_logistic_regression_rbd.R $1 6 161768456 163148994 result_rbd/PARK2_CNV.csv result_rbd/PARK2_CNV_control.csv

#Rscript build_logistic_regression_pd.R $1 1 155204073 155214819 result_v2/GBA_exon.csv result_v2/GBA_exon_control.csv

#Rscript build_logistic_regression_pd.R $1 4 90645125 90759594 result_v2/SNCA_exon.csv result_v2/SNCA_exon_control.csv

#Rscript build_logistic_regression_pd.R $1 1 20959900 20978175 result_v2/PINK1_exon.csv result_v2/PINK1_exon_control.csv

#Rscript build_logistic_regression_pd.R $1 12 40618681 40763272 result_v2/LRRK2_exon.csv result_v2/LRRK2_exon_control.csv

#Rscript build_logistic_regression_pd.R $1 16 46693467 46723323 result_v2/VPS35_exon.csv result_v2/VPS35_exon_control.csv
