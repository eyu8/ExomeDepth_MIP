#!/bin/bash

module load gcc/7.3.0 r/3.5.2
Rscript build_logistic_regression_pd.R $1 6 161768456 163148994 result_v$1/PARK2_exon_v$1.csv result_v$1/PARK2_exon_v$1_control.csv

Rscript build_logistic_regression_pd.R $1 1 155204073 155214819 result_v$1/GBA_exon_v$1.csv result_v$1/GBA_exon_v$1_control.csv

Rscript build_logistic_regression_pd.R $1 4 90645125 90759594 result_v$1/SNCA_exon_v$1.csv result_v$1/SNCA_exon_v$1_control.csv

Rscript build_logistic_regression_pd.R $1 1 20959900 20978175 result_v$1/PINK1_exon_v$1.csv result_v$1/PINK1_exon_v$1_control.csv

Rscript build_logistic_regression_pd.R $1 12 40618681 40763272 result_v$1/LRRK2_exon_v$1.csv result_v$1/LRRK2_exon_v$1_control.csv

Rscript build_logistic_regression_pd.R $1 16 46693467 46723323 result_v$1/VPS35_exon_v$1.csv result_v$1/VPS35_exon_v$1_control.csv
