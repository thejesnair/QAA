#!/bin/bash
#SBATCH --account=bgmp       
#SBATCH --partition=bgmp      
#SBATCH --job-name=QS_distr 
#SBATCH --output=QS_distr.out       #STDOUT
#SBATCH --error=QS_distr.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime
#SBATCH --nodes=1              ### Node count required for the job
#SBATCH --cpus-per-task=1    
#SBATCH --mem=16G               #increased from 16 so dont run out of memory


#Variables
Cco_1=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1.fastq.gz
Cco_2=/projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2.fastq.gz
CcoxCrh_1=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1.fastq.gz
CcoxCrh_2=/projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2.fastq.gz

Cco_1_png=Cco_com101_EO_adult_1_1_QSdistr.png
Cco_2_png=Cco_com101_EO_adult_1_2.QSdistr.png
CcoxCrh_1_png=CcoxCrh_comrhy110_EO_adult_1_1_QSdistr.png
CcoxCrh_2_png=CcoxCrh_comrhy110_EO_adult_1_2_QSdistr.png

PY_SCRIPT=/projects/bgmp/tnair/bioinfo/Bi623/PS2/QS_distr.py

/usr/bin/time -v $PY_SCRIPT -f $Cco_1 -o $Cco_1_png -l 150

/usr/bin/time -v $PY_SCRIPT -f $Cco_2 -o $Cco_2_png -l 150

/usr/bin/time -v $PY_SCRIPT -f $CcoxCrh_1 -o $CcoxCrh_1_png -l 150

/usr/bin/time -v $PY_SCRIPT -f $CcoxCrh_2 -o $CcoxCrh_2_png -l 150