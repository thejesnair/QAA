#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --job-name=STAR_alignment
#SBATCH --output=STAR_alignment.out       #STDOUT
#SBATCH --error=STAR_alignment.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime
#SBATCH --mem=16G  

#Variables

#Campie database directory

genome_dir=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.STAR_2.7.11b

#Cco
Cco_R1=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1_QUALTRIM_paired.fastq.gz
Cco_R2=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2_QUALTRIM_paired.fastq.gz
Cco_alignment=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_

#CcoxCrh
CcoxCrh_R1=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1_QUALTRIM_paired.fastq.gz
CcoxCrh_R2=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2_QUALTRIM_paired.fastq.gz
CcoxCrh_alignment=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_

conda activate QAA


/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
--outFilterMultimapNmax 3 \
--outSAMunmapped Within KeepPairs \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn  $Cco_R1 $Cco_R2 \
--genomeDir $genome_dir \
--outFileNamePrefix $Cco_alignment

/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
--outFilterMultimapNmax 3 \
--outSAMunmapped Within KeepPairs \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--readFilesCommand zcat \
--readFilesIn  $CcoxCrh_R1 $CcoxCrh_R2 \
--genomeDir $genome_dir \
--outFileNamePrefix $CcoxCrh_alignment