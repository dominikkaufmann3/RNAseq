export PATH=/data/users/lfalquet/SBC07107_24/scripts:/software/bin:/usr/local/bin:$PATH;
export PATH=/data/users/lfalquet/SBC07107_24/scripts/canu-2.2/bin:/data/users/lfalquet/SBC07104_23/scripts/Flye-2.9/bin:$PATH
export LANG=en_US.UTF-8

module add MultiQC/1.11-foss-2021a 
module load EMBOSS/6.6.0-foss-2021a
module load R-bundle-IBU/2023121400-foss-2021a-R-4.3.2
module load SeqKit/2.6.1
module load FastQC/0.11.9-Java-11
module load MultiQC/1.11-foss-2021a    
module load BWA/0.7.17-GCC-10.3.0 
module load SPAdes/3.15.3-GCC-10.3.0  
module load Bowtie2/2.4.4-GCC-10.3.0 
module load prokka/1.14.5-gompi-2021a
module load MAFFT/7.490-GCC-10.3.0-with-extensions
module load CD-HIT/4.8.1-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0 
module load GATK/4.2.6.1-GCCcore-10.3.0-Java-11 
module load fastp/0.23.4-GCC-10.3.0 
module load Trimmomatic/0.39-Java-11 
module load BLAST+/2.15.0-gompi-2021a
module load SAMtools/1.13-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0 
module load tabix/0.2.6-GCCcore-10.3.0
module load Roary/3.13.0-foss-2021a
module load fastp/0.23.4-GCC-10.3.0
export EMBOSS_FILTER=1



#module add UHTS/Analysis/vcftools/0.1.15;
#module add UHTS/Analysis/HTSeq/0.9.1;
#module add UHTS/Analysis/mummer/4.0.0beta1
#module add Phylogeny/FastTree/2.1.10
#module add UHTS/Analysis/picard-tools/2.21.8
#module add UHTS/Analysis/Roary/3.11.0
#module add UHTS/Analysis/minimap2/2.17
#module add UHTS/Analysis/busco/4.1.4
#module load Development/java/17.0.6
