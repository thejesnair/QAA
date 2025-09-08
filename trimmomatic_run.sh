#!/bin/bash
#SBATCH --account=bgmp       
#SBATCH --partition=bgmp      
#SBATCH --job-name=trimmomatic 
#SBATCH --output=trimmomatic.out       #STDOUT
#SBATCH --error=trimmomatic.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime
#SBATCH --nodes=1              ### Node count required for the job
#SBATCH --cpus-per-task=8    #multithreaded, so increased to 8
#SBATCH --mem=16G 

#Variables

#Input files
Cco_1_trimmed=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1_TRIMMED.fastq.gz
Cco_2_trimmed=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2_TRIMMED.fastq.gz

CcoxCrh_1_trimmed=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1_TRIMMED.fastq.gz
CcoxCrh_2_trimmed=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2_TRIMMED.fastq.gz

#Output files, paired

Cco_1_qualtrim_paired=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1_QUALTRIM_paired.fastq.gz
Cco_2_qualtrim_paired=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2_QUALTRIM_paired.fastq.gz

CcoxCrh_1_qualtrim_paired=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1_QUALTRIM_paired.fastq.gz
CcoxCrh_2_qualtrim_paired=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2_QUALTRIM_paired.fastq.gz

#Output files, unpaired
Cco_1_qualtrim_UNpaired=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1_QUALTRIM_UNpaired.fastq.gz
Cco_2_qualtrim_UNpaired=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2_QUALTRIM_UNpaired.fastq.gz

CcoxCrh_1_qualtrim_UNpaired=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1_QUALTRIM_UNpaired.fastq.gz
CcoxCrh_2_qualtrim_UNpaired=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2_QUALTRIM_UNpaired.fastq.gz

conda activate QAA

/usr/bin/time -v trimmomatic PE -threads 8 -phred33 $Cco_1_trimmed $Cco_2_trimmed \
    $Cco_1_qualtrim_paired $Cco_1_qualtrim_UNpaired \
    $Cco_2_qualtrim_paired $Cco_2_qualtrim_UNpaired \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35 

/usr/bin/time -v trimmomatic PE -threads 8 -phred33 $CcoxCrh_1_trimmed $CcoxCrh_2_trimmed \
    $CcoxCrh_1_qualtrim_paired $CcoxCrh_1_qualtrim_UNpaired \
    $CcoxCrh_2_qualtrim_paired $CcoxCrh_2_qualtrim_UNpaired \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:35