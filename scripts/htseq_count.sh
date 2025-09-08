#!/bin/bash
#SBATCH --account=bgmp       
#SBATCH --partition=bgmp      
#SBATCH --job-name=htseq_count2 
#SBATCH --output=htseq_count2.out       #STDOUT
#SBATCH --error=htseq_count2.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime
#SBATCH --nodes=1              ### Node count required for the job
#SBATCH --cpus-per-task=1   
#SBATCH --mem=16G  


#VARIABLES

Cco_bam=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.bam
campie_gtf=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.gtf
Cco_stranded_yes=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_counts_stranded_yes.txt
Cco_stranded_reverse=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_counts_stranded_reverse.txt

CcoxCrh_bam=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_dedup.bam
CcoxCrh_stranded_yes=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_counts_stranded_yes.txt
CcoxCrh_stranded_reverse=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_counts_stranded_reverse.txt



#Cco
/usr/bin/time -v htseq-count --format=bam --order=pos --stranded=yes --type=exon --idattr=gene_id $Cco_bam $campie_gtf > $Cco_stranded_yes
/usr/bin/time -v htseq-count --format=bam --order=pos --stranded=reverse --type=exon --idattr=gene_id $Cco_bam $campie_gtf > $Cco_stranded_reverse

#CcoxCrh
/usr/bin/time -v htseq-count --format=bam --order=pos --stranded=yes --type=exon --idattr=gene_id $CcoxCrh_bam $campie_gtf > $CcoxCrh_stranded_yes
/usr/bin/time -v htseq-count --format=bam --order=pos --stranded=reverse --type=exon --idattr=gene_id $CcoxCrh_bam $campie_gtf > $CcoxCrh_stranded_reverse