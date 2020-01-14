#!/bin/bash

parallel 'Rscript getCount_nobatch_Case.R {1} {2} {2} {3}' ::: $1 ::: $(seq $2 $3) ::: $4
