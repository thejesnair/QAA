#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --job-name=genome_database
#SBATCH --output=genome_database%j.out       #STDOUT
#SBATCH --error=genome_database%j.err       #STDERR
#SBATCH --time=0-4:00:00      ##Walltime
#SBATCH --mem=32G  


conda activate QAA #star in QAA env

/usr/bin/time -v STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.STAR_2.7.11b \
--genomeFastaFiles /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.fasta \
--sjdbGTFfile /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.gtf \
--genomeSAindexNbases 13