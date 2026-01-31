#!/bin/bash
#SBATCH --job-name=samtools_sort
#SBATCH --output=samtools_sort_%j.out
#SBATCH --error=samtools_sort_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=25000
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8

CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
MAP_DIR="/data/users/dkaufmann3/RNAseq/map"

for BAM in ${MAP_DIR}/*.bam
do
    SAMPLE=$(basename ${BAM} .bam)

    apptainer exec --bind /data ${CONTAINER} \
        samtools sort -@ ${SLURM_CPUS_PER_TASK} \
            -m 1G \
            -T ${MAP_DIR}/temp_${SAMPLE} \
            -o ${MAP_DIR}/${SAMPLE}.sorted.bam \
            ${BAM}

    # remove unsorted BAM to save space
    # rm ${BAM}
done
