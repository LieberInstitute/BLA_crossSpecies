#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=200G
#$ -pe local 8
#$ -N BLA_molprofile_cellranger
#$ -o logs/Human/cellranger_run1.$TASK_ID.txt
#$ -e logs/Human/cellranger_run1.$TASK_ID.txt
#$ -m e
#$ -t 1-5
#$ -tc 5

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load CellRanger
module load cellranger/7.0.0

## List current modules for reproducibility
module list

## Locate sample
SAMPLE=$(awk "NR==${SGE_TASK_ID}" 01_cellranger_human.txt)
echo "Processing sample ${SAMPLE}"
echo "${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=/dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/raw-data/FASTQ/Human/${SAMPLE} \
    --sample=${SAMPLE} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/processed-data/01_cellranger/Human/
mv ${SAMPLE} /dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/processed-data/01_cellranger/Human/

echo "**** Job ends ****"
date

