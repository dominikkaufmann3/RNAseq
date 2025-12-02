#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --output=hisat2_index_%j.out
#SBATCH --error=hisat2_index_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000M
#SBATCH --time=03:00:00
#SBATCH --partition=pibu_el8

CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
REF_DIR="/data/users/dkaufmann3/RNAseq/map"
INDEX_DIR="${REF_DIR}/genome_index"
INDEX_BASE="${INDEX_DIR}/genome"

mkdir -p ${REF_DIR} ${INDEX_DIR}

FA_URL="https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz"

FA_GZ="${REF_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
GTF_GZ="${REF_DIR}/Mus_musculus.GRCm39.115.gtf.gz"
FA="${REF_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF="${REF_DIR}/Mus_musculus.GRCm39.115.gtf"

# Download zipped files if not already present
if [[ ! -f ${FA_GZ} ]]; then
    wget -O ${FA_GZ} ${FA_URL}
fi

if [[ ! -f ${GTF_GZ} ]]; then
    wget -O ${GTF_GZ} ${GTF_URL}
fi

# Decompress for indexing if decompressed files don't exist
if [[ ! -f ${FA} ]]; then
    gunzip -c ${FA_GZ} > ${FA}
fi

if [[ ! -f ${GTF} ]]; then
    gunzip -c ${GTF_GZ} > ${GTF}
fi

# Build HISAT2 index
apptainer exec --bind /data ${CONTAINER} hisat2-build \
    -p ${SLURM_CPUS_PER_TASK} \
    ${FA} ${INDEX_BASE}
