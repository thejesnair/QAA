### PS2: RNA-seq Quality Assessment Assignment 

*GOAL:*<BR>
Process electric organ/skeletal muscle RNA-seq reads for eventual DGE analysis. Use tools for quality assessment, read trimming, compare quality assessments to those created by own software, and how to align and count reads.

Will also summarize important information in a high-level report. 

Make a cohesive report for your PI on what you've learned from data.

### Dataset

2 RNA-seq files, 2 diff electric fish studies (PRJNA1005245 and PRJNA1005244):<br>

-------------

- Download data from NCBI SRA
- Dump into FASTQ files
- Zip files (check ICA1 for refresh)

- Processing data for future assignment, keep files well organized.

### Rename files:

**Species_sample_tissue_age/size_sample#_readnumber.fastq.gz**<br><br>

Metadata (and filenames) can be found in "RNAseq_meta_campys.csv"<br>

SRR25630384	Thejes
- Rename: CcoxCrh_comrhy110_EO_adult_1_1.fastq.gz
- Rename: CcoxCrh_comrhy110_EO_adult_1_2.fastq.gz


SRR25630408	Thejes
- Rename: Cco_com101_EO_adult_1_1.fastq.gz
- Rename: Cco_com101_EO_adult_1_2.fastq.gz

```bash
mv SRR25630384_1.fastq.gz CcoxCrh_comrhy110_EO_adult_1_1.fastq.gz

```

# Part 1: Download data, Read QS Distributions

## First let's create conda env, QAA

Issues installing cutadapt, not compatible with current Python (anything older than 3.12)<br>
Had to install older version of Python (3.10) in QAA environment

```{bash}

srun -A bgmp -p bgmp --time=1:00:00 --pty bash
conda create -n QAA -c conda-forge -c bioconda python=3.10 fastqc=0.12.1 trimmomatic=0.39 cutadapt=5.0

fastqc --version
trimmomatic -version
cutadapt --version
```
FastQC: 0.12.1
Trimmomatic: 0.40
Cutadapt: 5.0

Initially installed trimmomatic 0.40, also needed to install older python so deleted QAA env and restarted
```{bash}
conda remove -n QAA --all
#deleted environment and then ran code above
```

Fetching SRA data

```{bash}
prefetch SRR25630408
prefetch SRR25630384
```

set up slurm script to do parallel runs for
fasterq-dump <br>
pigz: (leslie mentioned earlier over summer as better way to compress vs gzip, parallel implementation of gzip. already installed and available on Talapas)

SRA_fastq.sh
```bash
spots read      : 8,870,047
reads read      : 17,740,094
reads written   : 17,740,094
Command being timed: "fasterq-dump /projects/bgmp/tnair/bioinfo/Bi623/PS2/SRR25630384/SRR25630384.sra --threads 8"
User time (seconds): 35.44
System time (seconds): 5.06
Percent of CPU this job got: 427%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:09.48
Maximum resident set size (kbytes): 756692
Exit status: 0

spots read      : 54,900,823
reads read      : 109,801,646
reads written   : 109,801,646
Command being timed: "fasterq-dump /projects/bgmp/tnair/bioinfo/Bi623/PS2/SRR25630408/SRR25630408.sra --threads 8"
User time (seconds): 217.23
System time (seconds): 31.57
Percent of CPU this job got: 496%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:50.07
Maximum resident set size (kbytes): 756196
Exit status: 0

Command being timed: "pigz -p 8 SRR25630384_1.fastq SRR25630384_2.fastq"
User time (seconds): 561.11
System time (seconds): 0.81
Percent of CPU this job got: 747%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:15.15
Maximum resident set size (kbytes): 9560
Exit status: 0

Command being timed: "pigz -p 8 SRR25630408_1.fastq SRR25630408_2.fastq"
User time (seconds): 3455.61
System time (seconds): 5.39
Percent of CPU this job got: 749%
Elapsed (wall clock) time (h:mm:ss or m:ss): 7:41.54
Maximum resident set size (kbytes): 8704
Exit status: 0
```

## FastQC plots

- per base QS distributions for R1 and R2
    - "per_base_quality.png"
- per-base N content for R1 and R2
    - "per_base_n_content.png"
- consistent with QS plots?

#command line
```bash
conda activate QAA
#need to use fastqc which is installed in QAA env

fastqc SRR25630384_1.fastq.gz SRR25630384_2.fastq.gz
fastqc SRR25630408_1.fastq.gz SRR25630408_2.fastq.gz
zipinfo SRR25630384_1_fastqc.zip #view contents of zip file, know which figures to pull

unzip -p SRR25630384_1_fastqc.zip "*/Images/per_base_quality.png" > CcoxCrh_fastqc_per_base_quality_1.png #pull per base quality png, -p redirect files to stdout and save w/ new name, dont have to unzip entire .zip to get single file
unzip -p SRR25630384_1_fastqc.zip "*/Images/per_base_n_content.png" > CcoxCrh_fastqc_per_base_n_content_1.png #pull per base n content png
```

Note:
Could not do slurm for fastqc (maintenance, but learned later on have to pass --time flag in script: specify 0 days and however many hrs [did 4])
Decided to rerun SRR25630384 to get estimate of time for LN and to answer q's in assignment. 

fastq_rerun.sh slurm results

```bash
Command being timed: "fastqc /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1.fastq.gz /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2.fastq.gz -o /projects/bgmp/tnair/bioinfo/Bi623/PS2/fastq_rerun"
User time (seconds): 549.31
System time (seconds): 29.08
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 9:42.54
Maximum resident set size (kbytes): 593452
Exit status: 0
```

Run QS_distr.py script from Demux (Bi622)
Requires argparse arguments, file name, output png name, list length (length of sequence)
```bash
zcat Cco_com101_EO_adult_1_1.fastq.gz | head #determine length of sequences for QS_distr.py
#length = 150, stated in header!
```

QSdistr_usingdemuxscript.sh
```bash
Command being timed: "/projects/bgmp/tnair/bioinfo/Bi623/PS2/QS_distr.py -f /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1.fastq.gz -o Cco_com101_EO_adult_1_1_QSdistr_V2.png -l 150"
User time (seconds): 902.40
System time (seconds): 0.52
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 15:06.04
Maximum resident set size (kbytes): 66028
Exit status: 0

Command being timed: "/projects/bgmp/tnair/bioinfo/Bi623/PS2/QS_distr.py -f /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2.fastq.gz -o Cco_com101_EO_adult_1_2.QSdistr_V2.png -l 150"
User time (seconds): 880.46
System time (seconds): 0.63
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 14:48.96
Maximum resident set size (kbytes): 69720
Exit status: 0

Command being timed: "/projects/bgmp/tnair/bioinfo/Bi623/PS2/QS_distr.py -f /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1.fastq.gz -o CcoxCrh_comrhy110_EO_adult_1_1_QSdistr_V2.png -l 150"
User time (seconds): 141.71
System time (seconds): 0.13
Percent of CPU this job got: 97%
Elapsed (wall clock) time (h:mm:ss or m:ss): 2:24.83
Maximum resident set size (kbytes): 70452
Exit status: 0

Command being timed: "/projects/bgmp/tnair/bioinfo/Bi623/PS2/QS_distr.py -f /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2.fastq.gz -o CcoxCrh_comrhy110_EO_adult_1_2_QSdistr_V2.png -l 150"
User time (seconds): 141.43
System time (seconds): 0.13
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 2:21.99
Maximum resident set size (kbytes): 64204
Exit status: 0
```

Note: Figures have blocks in middle, essentially because bars in those positions have nearly identical heights/QS so is merging in the figure. No issue with data.

_V2, was purely to see if I could improve the look of the figure because of the darker blocks in the middle. Will most likely go with original figure. Not enough of a change to warrant adjusting script or going with version2.

Additional resource for fastqc assessment
https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html?utm_source=chatgpt.com

# Part 2: Adapter trimming comparison

Identified adapters looking at Adapter Content graph in .html file<br>
Additionally used command to look at Adapter content module in fastaqc_data.txt file<br>

Looking at the graph and data columns we can see that the Illumina Universal Adapter presence increases towards the end of the read,
whereas all the other possible adapters remain near 0. 

After further reseach, this adapter is the most common since most libraries use TruSeq adapters.

Read documentation on adapters.
https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm

Adapters located towards 3' end.

R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

```bash
unzip -p SRR25630384_1_fastqc.zip "*/fastqc_data.txt" | awk '/>>Adapter Content/,/>>END_MODULE/' | sed '1d;$d'

#Determine occurrences of adapter in R1,R2, checking for adapters in 3' tail
zcat CcoxCrh_comrhy110_EO_adult_1_1.fastq.gz | sed -n '2~4p' | awk '{print substr($0,length($0)-29)}' | grep -Ec 'AGATCGGAAGAGC'
zcat CcoxCrh_comrhy110_EO_adult_1_2.fastq.gz | sed -n '2~4p' | awk '{print substr($0,length($0)-29)}' | grep -Ec 'AGATCGGAAGAGC'

#Check adapter/sequence orientation, checking for adapter seq in 5' start (should be very low)
zcat Cco_com101_EO_adult_1_1.fastq.gz | sed -n '2~4p' | awk '{print substr($0,1,30)}' | grep -Ec '^AGATCGGAAGAGC'
zcat Cco_com101_EO_adult_1_2.fastq.gz | sed -n '2~4p' | awk '{print substr($0,1,30)}' | grep -Ec '^AGATCGGAAGAGC'
```
3' end
Cco_1: 2533728
Cco_2: 2568303

CcoxCrh_1: 627873
CcoxCrh_2: 624953

5' end
Cco_1: 1
Cco_2: 1

CcoxCrh_1: 0
CcoxCrh_2: 0

### Cutadapt
https://cutadapt.readthedocs.io/en/v5.0/reference.html#command-line-options

Slurm script to run cutadapt on fastq.gz files<br>
Using logs will output stdout to a file for r1 and r2, and stderr to another file. 
``` > logs/cutadapt_CcoxCrh.txt 2>&1 ```


cutadapt_run.sh
```bash
Cco_com101_EO_adult_1
  R1: 17.4%
  R2: 17.0%

CcoxCrh_comrhy110_EO_adult_1
  R1: 21.8%
  R2: 21.6%
```

```
Command being timed: "cutadapt -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1_TRIMMED.fastq.gz -p /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2_TRIMMED.fastq.gz /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1.fastq.gz /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2.fastq.gz"
User time (seconds): 482.12
System time (seconds): 54.95
Percent of CPU this job got: 478%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:52.18
Maximum resident set size (kbytes): 93640
Exit status: 0

Command being timed: "cutadapt -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1_TRIMMED.fastq.gz -p /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2_TRIMMED.fastq.gz /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1.fastq.gz /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2.fastq.gz"
User time (seconds): 85.03
System time (seconds): 8.22
Percent of CPU this job got: 501%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:18.59
Maximum resident set size (kbytes): 85384
Exit status: 0
```

### Trimmomatic

Run trimmomatic on the cutadapt trimmed reads, we are quality trimming the reads now. (Though trimmomatic is capable of adapter trimming as well)

https://gensoft.pasteur.fr/docs/Trimmomatic/0.39/TrimmomaticManual_V0.32.pdf


Initial look into some of the different read length distributions
```bash
zcat Cco_com101_EO_adult_1_1_QUALTRIM_paired.fastq.gz | awk 'NR%4==2 { print length($0) }' | head -n 20
```

plot_trimmed_distr.py and plot_trimmed_distr.sh

```bash
Command being timed: "./plot_trimmed_distr.py -R1 /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1_QUALTRIM_paired.fastq.gz -R2 /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2_QUALTRIM_paired.fastq.gz -o Cco_read_len_distr.png -t Cco_com101_EO_adult_1 Trimmed Read Length Distribution"
User time (seconds): 181.69
System time (seconds): 0.70
Percent of CPU this job got: 98%
Elapsed (wall clock) time (h:mm:ss or m:ss): 3:05.89
Maximum resident set size (kbytes): 907108
Exit status: 0

Command being timed: "./plot_trimmed_distr.py -R1 /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1_QUALTRIM_paired.fastq.gz -R2 /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2_QUALTRIM_paired.fastq.gz -o CcoxCrh_read_len_distr.png -t CcoxCrh_comrhy101_EO_adult_1 Trimmed Read Length Distribution"
User time (seconds): 30.72
System time (seconds): 0.21
Percent of CPU this job got: 98%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:31.27
Maximum resident set size (kbytes): 204624
Exit status 0
```
# Part 3: Alignment and strand-specificity

Bioconda install:
- Star
- Picard
- Samtools
- NumPy
- Matplotlib
- HTSeq

```bash
conda activate QAA

conda install bioconda::star
conda install bioconda::picard
conda install bioconda::samtools
conda install bioconda::htseq

conda install -c conda-forge numpy matplotlib

#verify installation
STAR --version
samtools --version
conda list picard
htseq-count --version
python -c  "import numpy, matplotlib; print(numpy.__version__); print(matplotlib.__version__)"
```

STAR: 2.7.11b<br>
samtools:1.22.1<br>
picard 3.4.0<br>
htseq-count 2.0.9<br>
numpy 1.26.4<br>
matplotlib 3.10.6

Download Campie data from Dryad<br>
Unable to use wget, Error 403, forbidden. Clicking on link works but can't download to Talapas<br>
Copying files from Talapas location instead

```bash
cp /projects/bgmp/shared/Bi623/PS2/campylomormyrus.fasta .
cp /projects/bgmp/shared/Bi623/PS2/campylomormyrus.gff .

#convert gff to gtf
conda install bioconda::gffread
gffread campylomormyrus.gff -T -o campylomormyrus.gtf
```

Used PS8_genome.srun script to generate STAR database<br>
Changed name to genome_database_PS8script.sh


Ran script and got error

```bash
!!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=862592683, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 13
/gpfs/home/tnair/miniforge3/envs/QAA/bin/STAR: line 8: 3987328 Killed                  "${cmd}" "$@"
Command exited with non-zero status 137
```

Memory error based on default setting of --genomeSAindexbases. Lowered parameter to 13 based on suggestion and included flag in .sh<br>
Another error encountered: essentially ran out of memory<br>
Increased mem to 32, STAR index requires more memory? (Even though Campie assembly is smaller, it may not be as good quality as Ensembl's zebrafish, might have
more contigs/scaffolds/repeats so STAR eating up more RAM trying to build a more complex index)


```bash
slurmstepd: error: Detected 1 oom_kill event in StepId=38058703.batch. Some of the step tasks have been OOM Killed.
```
- increasing memory to 32 worked

genome_database38058952.err
```bash
Command being timed: "STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.STAR_2.7.11b --genomeFastaFiles /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.fasta --sjdbGTFfile /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.gtf --genomeSAindexNbases 13"
User time (seconds): 1518.07
System time (seconds): 23.36
Percent of CPU this job got: 453%
Elapsed (wall clock) time (h:mm:ss or m:ss): 5:39.69
Maximum resident set size (kbytes): 26030040
Exit status 0 
```

Align reads to STAR campie database. Used same settings from STAR_alignment.sh from PS8. 

STAR_alignment.sh
```bash
Command being timed: "STAR --runThreadN 8 --runMode alignReads --outFilterMultimapNmax 3 --outSAMunmapped Within KeepPairs --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesIn /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_1_QUALTRIM_paired.fastq.gz /projects/bgmp/tnair/bioinfo/Bi623/PS2/Cco_com101_EO_adult_1_2_QUALTRIM_paired.fastq.gz --genomeDir /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.STAR_2.7.11b --outFileNamePrefix /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_"
User time (seconds): 9180.89
System time (seconds): 20.26
Percent of CPU this job got: 773%
Elapsed (wall clock) time (h:mm:ss or m:ss): 19:49.44
Maximum resident set size (kbytes): 10256084
Exit status 0

Command being timed: "STAR --runThreadN 8 --runMode alignReads --outFilterMultimapNmax 3 --outSAMunmapped Within KeepPairs --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesIn /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_1_QUALTRIM_paired.fastq.gz /projects/bgmp/tnair/bioinfo/Bi623/PS2/CcoxCrh_comrhy110_EO_adult_1_2_QUALTRIM_paired.fastq.gz --genomeDir /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.STAR_2.7.11b --outFileNamePrefix /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_"
User time (seconds): 1291.21
System time (seconds): 6.49
Percent of CPU this job got: 739%
Elapsed (wall clock) time (h:mm:ss or m:ss): 2:55.60
Maximum resident set size (kbytes): 10076004
Exit status 0
```

Convert to BAM and remove PCR duplication with Picard
https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

deduplication.sh 

Error because didn't pass read groups<br> 
Sam conversion exited with 0 issues<br>
Need to rerun for read groups and picard, will comment out sam conversion in slurm script.<br>

```bash
Command being timed: "samtools view -b /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_Aligned.out.sam"
User time (seconds): 479.41
System time (seconds): 5.72
Percent of CPU this job got: 64%
Elapsed (wall clock) time (h:mm:ss or m:ss): 12:28.77
Maximum resident set size (kbytes): 3648
Exit status 0

Command being timed: "samtools view -b /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_Aligned.out.sam"
User time (seconds): 83.69
System time (seconds): 0.87
Percent of CPU this job got: 76%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:50.22
Maximum resident set size (kbytes): 3700
Exit status 0

Command exited with non-zero status 1
Command being timed: "picard MarkDuplicates INPUT=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_sorted.bam OUTPUT=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.bam METRICS_FILE=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT"
User time (seconds): 1.30
System time (seconds): 0.18
Percent of CPU this job got: 53%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.76
Maximum resident set size (kbytes): 165728
Exit status: 1

Command exited with non-zero status 1
Command being timed: "picard MarkDuplicates INPUT=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_sorted.bam OUTPUT=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_dedup.bam METRICS_FILE=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_dedup.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT"
User time (seconds): 1.32
System time (seconds): 0.13
Percent of CPU this job got: 98%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.47
Maximum resident set size (kbytes): 168724
Exit status: 1
```

rerun with read group addition
```bash
Command being timed: "picard AddOrReplaceReadGroups I=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_sorted.bam O=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_sorted_RG.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=Cco VALIDATION_STRINGENCY=LENIENT"
User time (seconds): 588.86
System time (seconds): 4.82
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 9:57.67
Maximum resident set size (kbytes): 265380
Exit status 0

Command being timed: "picard AddOrReplaceReadGroups I=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_sorted.bam O=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_sorted_RG.bam RGID=2 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=CcoxCrh VALIDATION_STRINGENCY=LENIENT"
User time (seconds): 107.04
System time (seconds): 1.06
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:48.81
Maximum resident set size (kbytes): 265032
Exit status: 0

Command being timed: "picard MarkDuplicates INPUT=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_sorted_RG.bam OUTPUT=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.bam METRICS_FILE=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT"
User time (seconds): 697.28
System time (seconds): 9.32
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 11:51.12
Maximum resident set size (kbytes): 2188316
Exit status 0

Command being timed: "picard MarkDuplicates INPUT=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_sorted_RG.bam OUTPUT=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_dedup.bam METRICS_FILE=/projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_dedup.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT"
User time (seconds): 135.32
System time (seconds): 2.35
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 2:18.61
Maximum resident set size (kbytes): 2159036
Exit status: 0
```

Need to run SAM files using is_mapped.py script from PS8 to determine mapped and unmapped read counts.<br>
Files are currently in BAM format. 

Commands ran in terminal
```bash
#convert to SAM
samtools view -h -o Cco_dedup.sam Cco_dedup.bam
samtools view -h -o CcoxCrh_dedup.sam CcoxCrh_dedup.bam

#Count reads
./is_mapped.py -f /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.sam
```

HTSEQ count deduplicated reads
https://htseq.readthedocs.io/en/release_0.11.1/count.html
 
```bash
#check what the gene id col looks like for --idattr (-i) parameter
head campylomormyrus.gtf

#Cco
/usr/bin/time -v htseq-count --format=bam --order=pos --stranded=yes --type=exon --idattr=gene_id $Cco_bam $campie_gtf > $Cco_stranded_yes
/usr/bin/time -v htseq-count --format=bam --order=pos --stranded=reverse --type=exon --idattr=gene_id $Cco_bam $campie_gtf > $Cco_stranded_reverse
```
htseq_count.sh
```bash
Command being timed: "htseq-count --format=bam --order=pos --stranded=yes --type=exon --idattr=gene_id /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.bam /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.gtf"
User time (seconds): 1212.93
System time (seconds): 2.08
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 20:23.25
Maximum resident set size (kbytes): 2321684
Exit status: 0

Command being timed: "htseq-count --format=bam --order=pos --stranded=reverse --type=exon --idattr=gene_id /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/Cco_dedup.bam /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.gtf"
User time (seconds): 1242.40
System time (seconds): 2.11
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 20:53.77
Maximum resident set size (kbytes): 2315564
Exit status: 0

Command being timed: "htseq-count --format=bam --order=pos --stranded=yes --type=exon --idattr=gene_id /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_dedup.bam /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.gtf"
User time (seconds): 278.83
System time (seconds): 0.43
Percent of CPU this job got: 98%
Elapsed (wall clock) time (h:mm:ss or m:ss): 4:44.67
Maximum resident set size (kbytes): 188208
Exit status: 0

Command being timed: "htseq-count --format=bam --order=pos --stranded=reverse --type=exon --idattr=gene_id /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/CcoxCrh_dedup.bam /projects/bgmp/tnair/bioinfo/Bi623/PS2/part3_alignment/campylomormyrus.gtf"
User time (seconds): 288.25
System time (seconds): 0.40
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 4:50.
Exit status: 0
```
Investigate the counts for strandedness
```bash
#determine no of assigned reads
awk '$1 !~ /^__/ {s+=$2} END{print s+0}' Cco_counts_stranded_yes.txt

#determine no feature, ambiguous reads 
tail Cco_counts_stranded_yes.txt
```

Cco_stranded_yes
assigned reads: 555332
__no_feature    16650247
__ambiguous     2263
__too_low_aQual 119120
__not_aligned   8853141
__alignment_not_unique  593473

Cco_stranded_reverse
assigned reads: 10328955
__no_feature    6747807
__ambiguous     131080
__too_low_aQual 119120
__not_aligned   8853141
__alignment_not_unique  593473

CcoxCrh_stranded_yes
assigned reads: 138529
__no_feature    4448568
__ambiguous     427
__too_low_aQual 3126
__not_aligned   491389
__alignment_not_unique  173254

CcoxCrh_stranded_reverse
assigned reads: 2745246
__no_feature    1808597
__ambiguous     33681
__too_low_aQual 3126
__not_aligned   491389
__alignment_not_unique  173254










### HTSeq-count results: stranded=yes vs stranded=reverse

| Sample    | Strandedness | Assigned Reads | __no_feature | __ambiguous |
|-----------|--------------|----------------|--------------|-------------|
| **Cco**   | yes          |     555,332    | 16,650,247   |     2,263   |
|           | reverse      |  10,328,955    |  6,747,807   |   131,080   |
| **CcoxCrh** | yes        |     138,529    |  4,448,568   |       427   |
|           | reverse      |   2,745,246    |  1,808,597   |    33,681   |