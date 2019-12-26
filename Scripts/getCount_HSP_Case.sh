#!/bin/bash

parallel 'Rscript getCount_HSP_Case.R {} {}' ::: $(seq $1 $2)
