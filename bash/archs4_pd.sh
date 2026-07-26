#!/bin/bash
#SBATCH --job-name=archs4_pd
#SBATCH --partition=cpu
#SBATCH --output=logs/archs4_pd_%A_%a.out
#SBATCH --error=logs/archs4_pd_%A_%a.err
#SBATCH --array=5
#SBATCH --cpus-per-task=6
#SBATCH --mem=90G
#SBATCH --time=22:00:00

mkdir -p logs

R_SCRIPT="/home/ctools/opt/R-4.4.1/bin/Rscript"

DATASETS=(AD CVD OB T2D COPD IBD RA TB)
DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"

echo "Running dataset: ${DATASET}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"


$R_SCRIPT ../R/01_run_archs4_pairedGSEA.R ${DATASET}
