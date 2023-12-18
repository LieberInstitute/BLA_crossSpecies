#!/bin/bash
#SBATCH --job-name=Norm_and_DimRed
#SBATCH --output=logs/Norm_and_DimRed.txt
#SBATCH --error=logs/Norm_and_DimRed.txt
#SBATCH --mem=200G
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=20 # specify the number of CPUs needed for the job, adjust as needed

echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


## Load the R module
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript norm_and_dim_reduction_all_primates.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/