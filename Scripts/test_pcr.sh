#!/bin/bash

Rscript test_pcr.R result_v$1/PARK2_exon_v$1.csv original.targets_only.sorted.bed

Rscript test_pcr.R result_v$1/GBA_exon_v$1.csv original.targets_only.sorted.bed

Rscript test_pcr.R result_v$1/SNCA_exon_v$1.csv original.targets_only.sorted.bed

Rscript test_pcr.R result_v$1/PINK1_exon_v$1.csv original.targets_only.sorted.bed

Rscript test_pcr.R result_v$1/LRRK2_exon_v$1.csv original.targets_only.sorted.bed

Rscript test_pcr.R result_v$1/VPS35_exon_v$1.csv original.targets_only.sorted.bed
