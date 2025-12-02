#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="bwa"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=36:00:00
#SBATCH --mem=25G
#SBATCH --partition=pcoursea

XX=$1
GROUP=$2
THREADS=$SLURM_CPUS_PER_TASK

#create and go to the TP directory
mkdir bwa$1
cd bwa$1

source /data/users/lfalquet/SBC07107_25/scripts/module.sh

#link the reference genome and the reads locally
ln -s /data/users/lfalquet/SBC07107_25/YEAST/reference/R64-1-1.112.fa .
ln -s /data/users/lfalquet/SBC07107_25/YEAST/reference/R64-1-1.112.gtf .
ln -s /data/users/lfalquet/SBC07107_25/YEAST/${GROUP}/${XX}_1.fastq.gz ${XX}_1.fastq.gz
ln -s /data/users/lfalquet/SBC07107_25/YEAST/${GROUP}/${XX}_2.fastq.gz ${XX}_2.fastq.gz

#ln -s /data/users/lfalquet/SBC07107_25/preparation/rawreads/${XX}_1.fq.gz ${XX}_1.fastq.gz
#ln -s /data/users/lfalquet/SBC07107_25/preparation/rawreads/${XX}_2.fq.gz ${XX}_2.fastq.gz

#verify that your files are in your directory (ls -l)

fastqc -t $THREADS ${XX}_1.fastq.gz ${XX}_2.fastq.gz

#with trimmotatic
#trimmomatic PE -threads ${THREADS} -phred33 ${XX}_1.fastq.gz ${XX}_2.fastq.gz ${XX}_1trim.fastq.gz ${XX}_1unpaired.fastq.gz ${XX}_2trim.fastq.gz ${XX}_2unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:130

#with fastp
fastp -i ${XX}_1.fastq.gz -I ${XX}_2.fastq.gz -o ${XX}_1trim.fastq.gz -O ${XX}_2trim.fastq.gz -j ${XX}_fastp.json -h ${XX}_fastp.html --thread ${THREADS} --trim_poly_g -l 50;


#index reference for bwa
bwa index R64-1-1.112.fa

#map reads onto reference
bwa mem -t $THREADS -M R64-1-1.112.fa ${XX}_1trim.fastq.gz ${XX}_2trim.fastq.gz > ${XX}.sam

#index reference for samtools
samtools faidx R64-1-1.112.fa

#convert sam to bam
samtools view -b -@ $THREADS -t R64-1-1.112.fa.fai ${XX}.sam > ${XX}_unsorted.bam

#sort bam
samtools sort -@ $THREADS -o ${XX}.bam ${XX}_unsorted.bam
rm ${XX}_unsorted.bam  
rm ${XX}.sam

#index bam
samtools index -@ $THREADS ${XX}.bam

