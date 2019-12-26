#!/bin/bash

parallel 'Rscript getCount_HSP_non_affected_fam.R {} {}' ::: $(seq $1 $2)
