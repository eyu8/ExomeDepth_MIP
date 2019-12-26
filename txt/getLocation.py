#!/usr/bin/python

import sys

with open("MLPA_filtered.raw", "r") as raw:
    with open("MLPA_filtered.txt","w") as output:
        for line1 in raw:
            line = line1.split('\n',1)[0]
            path, number, name = line.split(':', 3)
            if path == "/home/eyu8/runs/eyu8/data/MIP/pd_noGBAP1_mip_filtered/pd.PD_MIP15.ibam.noGBAP1.name":
                cohort = 1
            elif path == "/home/eyu8/runs/eyu8/data/MIP/pd_noGBAP1_mip_filtered/pd.PD_MIP19.bam.noGBAP1.name":
                cohort = 2
            elif path == "/home/eyu8/runs/eyu8/data/MIP/pd_noGBAP1_mip_filtered/pd.PD_MIP22.bam.noGBAP1.name":
                cohort = 3
            elif path == "/home/eyu8/runs/eyu8/data/MIP/pd_noGBAP1_mip_filtered/pd.PD_MIP26.bam.noGBAP1.name":
                cohort = 4
            elif path == "/home/eyu8/runs/eyu8/data/MIP/pd_noGBAP1_mip_filtered/pd.PD_MIP28.bam.noGBAP1.name":
                cohort = 5
            elif path == "/home/eyu8/runs/eyu8/data/MIP/pd_noGBAP1_mip_filtered/pd.PD_MIP55.bam.noGBAP1.name":
                cohort = 6
            elif path == "/home/eyu8/runs/eyu8/data/MIP/pd_noGBAP1_mip_filtered/pd.PD_MIP65.bam.noGBAP1.name":
                cohort = 7
            elif path == "/home/eyu8/runs/eyu8/data/MIP/pd_noGBAP1_mip_filtered/pd.PD_MIP71.bam.noGBAP1.name":
                cohort = 8
            elif path == "/home/eyu8/runs/eyu8/data/MIP/pd_noGBAP1_mip_filtered/pd.PD_MIP76.bam.noGBAP1.name":
                cohort = 9
            output.write(name+' '+str(cohort)+' '+number+'\n')
