#!/bin/bash
#SBATCH --job-name=archs4_pd
#SBATCH --partition=gpu
#SBATCH --output=logs/archs4_pd_%A_%a.out
#SBATCH --error=logs/archs4_pd_%A_%a.err
#SBATCH --array=0-7
#SBATCH --cpus-per-task=8
#SBATCH --mem=45G
#SBATCH --time=48:00:00

mkdir -p logs

R_SCRIPT="/home/ctools/opt/R-4.4.1/bin/Rscript"

DATASETS=(AD CVD OB T2D COPD IBD RA TB)
DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"

echo "Running dataset: ${DATASET}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"


$R_SCRIPT ../R/01_run_archs4_pairedGSEA.R ${DATASET}
