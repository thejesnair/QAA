#!/bin/bash
#SBATCH --account=bgmp       
#SBATCH --partition=bgmp      
#SBATCH --job-name=cutadapt 
#SBATCH --output=cutadapt.out       #STDOUT
#SBATCH --error=cutadapt.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime
#SBATCH --nodes=1              ### Node count required for the job
#SBATCH --cpus-per-task=8    #multithreaded, so increased to 8
#SBATCH --mem=16G               


#Variables

#adapters
R1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
R2adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

#input files
Cco_1=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1.fastq.gz
Cco_2=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2.fastq.gz

CcoxCrh_1=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1.fastq.gz
CcoxCrh_2=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2.fastq.gz

#output files
Cco_1_output=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1_TRIMMED.fastq.gz
Cco_2_output=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2_TRIMMED.fastq.gz

CcoxCrh_1_output=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1_TRIMMED.fastq.gz
CcoxCrh_2_output=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2_TRIMMED.fastq.gz

conda activate QAA

echo "Running Cco trimming"
/usr/bin/time -v cutadapt -j 8 -a $R1adapter -A $R2adapter -o $Cco_1_output -p $Cco_2_output $Cco_1 $Cco_2

echo "Running CcoxCrh trimming"
/usr/bin/time -v cutadapt -j 8 -a $R1adapter -A $R2adapter -o $CcoxCrh_1_output -p $CcoxCrh_2_output $CcoxCrh_1 $CcoxCrh_2
