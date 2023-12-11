#!/bin/bash

#SBATCH --job-name=Baboon_cellranger
#SBATCH --output=logs/Baboon/cellranger_run1_ncbi.%A_%a.txt
#SBATCH --error=logs/Baboon/cellranger_run1_ncbi.%A_%a.txt
#SBATCH --cpus-per-task=8
#SBATCH --mem=70G
#SBATCH --array=1-7%7


echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
## load CellRanger
module load cellranger

## List current modules for reproducibility
module list

## Locate sample
SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" 01_cellranger_baboon.txt)
echo "Processing sample ${SAMPLE}"
echo "${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/raw-data/refdata/Panubis1_genome_ncbi \
    --fastqs=/dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/raw-data/FASTQ/Baboon/${SAMPLE} \
    --sample=${SAMPLE} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/processed-data/01_cellranger/Baboon_ncbi/
mv ${SAMPLE} /dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/processed-data/01_cellranger/Baboon_ncbi/

echo "**** Job ends ****"
date