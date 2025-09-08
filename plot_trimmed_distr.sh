#!/bin/bash
#SBATCH --account=bgmp       
#SBATCH --partition=bgmp      
#SBATCH --job-name=plot_readlen 
#SBATCH --output=plot_readlen.out       #STDOUT
#SBATCH --error=plot_readlen.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime
#SBATCH --nodes=1              ### Node count required for the job
#SBATCH --cpus-per-task=1   
#SBATCH --mem=16G  

#Variables
#Cco
Cco_R1=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1_QUALTRIM_paired.fastq.gz
Cco_R2=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2_QUALTRIM_paired.fastq.gz
Cco_out=Cco_read_len_distr.png
Cco_title="Cco_com101_EO_adult_1 Trimmed Read Length Distribution"

#CcoxCrh
CcoxCrh_R1=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1_QUALTRIM_paired.fastq.gz
CcoxCrh_R2=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2_QUALTRIM_paired.fastq.gz
CcoxCrh_out=CcoxCrh_read_len_distr.png
CcoxCrh_title="CcoxCrh_comrhy101_EO_adult_1 Trimmed Read Length Distribution"


/usr/bin/time -v ./plot_trimmed_distr.py -R1 $Cco_R1 -R2 $Cco_R2 -o $Cco_out -t "$Cco_title" 

/usr/bin/time -v ./plot_trimmed_distr.py -R1 $CcoxCrh_R1 -R2 $CcoxCrh_R2 -o $CcoxCrh_out -t "$CcoxCrh_title" 


