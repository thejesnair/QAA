#!/bin/bash
#SBATCH --account=bgmp       
#SBATCH --partition=bgmp      
#SBATCH --job-name=SRA_fastq 
#SBATCH --output=SRA_fastq.out       #STDOUT
#SBATCH --error=SRA_fastq.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime
#SBATCH --nodes=1              ### Node count required for the job
#SBATCH --cpus-per-task=8     #multithreaded, so increased to 8
#SBATCH --mem=32G               #increased from 16 so dont run out of memory

#activate conda environment
conda activate sra

#variables for SRR files
SRfile1=/projects/bgmp/tnair/bioinfo/Bi623/PS2/SRR25630384/SRR25630384.sra
SRfile2=/projects/bgmp/tnair/bioinfo/Bi623/PS2/SRR25630408/SRR25630408.sra

#fasterq-dump
echo "Staring fasterq-dump SRfile1"
/usr/bin/time -v fasterq-dump $SRfile1 --threads 8

echo "Starting fasterq-dump for SRfile2"
/usr/bin/time -v fasterq-dump $SRfile2 --threads 8 

#zip files
echo "Copressing FASTQ files"
/usr/bin/time -v pigz -p 8 SRR25630384_*.fastq
/usr/bin/time -v pigz -p 8 SRR25630408_*.fastq