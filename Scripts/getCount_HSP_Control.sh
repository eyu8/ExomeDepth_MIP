#!/bin/bash

parallel 'Rscript getCount_HSP_Control.R {} {}' ::: $(seq $1 $2)
