#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="snps"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=28:00:00
#SBATCH --mem=20G
#SBATCH --partition=pcoursea

NAME=$1
THREADS=$SLURM_CPUS_PER_TASK

source /data/users/lfalquet/SBC07107_25/scripts/module.sh

mkdir /data/users/${USER}/snp${NAME}
cd /data/users/${USER}/snp${NAME}

ln -s /data/users/lfalquet/SBC07107_25/YEAST/reference/R64-1-1.112.fa R64-1-1.112.fa
ln -s /data/users/lfalquet/SBC07107_25/YEAST/reference/R64-1-1.112.gtf R64-1-1.112.gtf
ln -s /data/users/${USER}/bwa*/*.bam* .

for SAMPLE in `ls *.bam`
do
bcftools mpileup -Ou --threads $THREADS -f R64-1-1.112.fa ${SAMPLE} | bcftools call -vc -Oz --threads $THREADS -o ${SAMPLE}.vcf.gz
tabix ${SAMPLE}.vcf.gz
done

bcftools merge -O v `ls *.vcf.gz` > ${NAME}all.vcf

SNPEFF=/data/users/lfalquet/SBC07107_25/scripts/snpEff
java -Xmx20g -jar ${SNPEFF}/snpEff.jar R64-1-1.105 -no-upstream -no-downstream ${NAME}all.vcf > ${NAME}all_annot.vcf 

cat ${NAME}all_annot.vcf | java -jar ${SNPEFF}/SnpSift.jar filter "(( ANN[*].EFFECT != 'synonymous_variant') & ( ANN[*].EFFECT != 'intergenic_region'))"   > ${NAME}all_coding.vcf

# keep only variants that are found in less than 3 strains
cat ${NAME}all_coding.vcf | java -jar ${SNPEFF}/SnpSift.jar filter  "((countVariant() < 3))" > ${NAME}all_coding_max2.vcf 

# Tip: modify before importing into excel
perl -pe 's/;ANN=/;\tANN=/g' ${NAME}all_coding_max2.vcf > ${NAME}all_coding_max2_forexcel.vcf 
perl -i -pe 's/FORMAT\t/FORMAT\t\t/' ${NAME}all_coding_max2_forexcel.vcf

