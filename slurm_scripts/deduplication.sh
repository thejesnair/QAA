#!/bin/bash
#SBATCH --account=bgmp       
#SBATCH --partition=bgmp      
#SBATCH --job-name=deduplication2 
#SBATCH --output=deduplication2.out       #STDOUT
#SBATCH --error=deduplication2.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime
#SBATCH --nodes=1              ### Node count required for the job
#SBATCH --cpus-per-task=1   
#SBATCH --mem=16G  


conda activate QAA

#VARIABLES
#SAM files
Cco_sam=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_Aligned.out.sam
CcoxCrh_sam=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_Aligned.out.sam

#SAM -> BAM files
Cco_sorted_bam=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_sorted.bam
CcoxCrh_sorted_bam=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_sorted.bam

#read groups added
Cco_sorted_RG=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_sorted_RG.bam
CcoxCrh_sorted_RG=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_sorted_RG.bam

#deduplication files, picard
Cco_dedup_bam=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.bam
CcoxCrh_dedup_bam=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_dedup.bam

#deduplication metrics, picard
Cco_dedup_metrics=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.metrics
CcoxCrh_dedup_metrics=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_dedup.metrics

#samtools convert SAM -> BAM
/usr/bin/time -v samtools view -b $Cco_sam | samtools sort -o $Cco_sorted_bam
/usr/bin/time -v samtools view -b $CcoxCrh_sam | samtools sort -o $CcoxCrh_sorted_bam

#Add read groups, because error in executing picard later on
/usr/bin/time -v picard AddOrReplaceReadGroups \
    I=$Cco_sorted_bam \
    O=$Cco_sorted_RG \
    RGID=1 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=Cco \
    VALIDATION_STRINGENCY=LENIENT

/usr/bin/time -v picard AddOrReplaceReadGroups \
    I=$CcoxCrh_sorted_bam \
    O=$CcoxCrh_sorted_RG \
    RGID=2 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=CcoxCrh \
    VALIDATION_STRINGENCY=LENIENT

#Deduplication with Picard
/usr/bin/time -v picard MarkDuplicates INPUT=$Cco_sorted_RG OUTPUT=$Cco_dedup_bam METRICS_FILE=$Cco_dedup_metrics \
REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT

/usr/bin/time -v picard MarkDuplicates INPUT=$CcoxCrh_sorted_RG OUTPUT=$CcoxCrh_dedup_bam METRICS_FILE=$CcoxCrh_dedup_metrics \
REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT