#!/bin/bash
#SBATCH --job-name=Baboon_build_sce
#SBATCH --output=logs/build_sce_baboon.txt
#SBATCH --error=logs/build_sce_baboon.txt
#SBATCH --mem=80G
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=mtotty2@jh.edu

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
Rscript build_sce_baboon.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/