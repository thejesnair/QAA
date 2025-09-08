#!/bin/bash
#SBATCH --account=bgmp       
#SBATCH --partition=bgmp      
#SBATCH --job-name=fastqc_rerun
#SBATCH --output=fastqc_rerun.out       #STDOUT
#SBATCH --error=fastqc_rerun.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime

conda activate QAA
/usr/bin/time -v fastqc /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1.fastq.gz /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2.fastq.gz -o /projects/bgmp/tnair/bioinfo/Bi623/PS2/fastq_rerun