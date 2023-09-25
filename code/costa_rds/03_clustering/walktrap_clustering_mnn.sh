#!/bin/bash
#$ -cwd
#$ -l mem_free=128G,h_vmem=128G,h_fsize=50G
#$ -pe local 4
#$ -N dimred_and_clustering
#$ -o logs/walktrap_clustering_mnn.txt
#$ -e logs/walktrap_clustering_mnn.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript walktrap_clustering_mnn.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/

