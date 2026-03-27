#!/bin/bash
#SBATCH --job-name=tcga_BRCA_LumA
#SBATCH --partition=gpu
#SBATCH --output=logs/tcga_BRCA_LumA_%j.out
#SBATCH --error=logs/tcga_BRCA_LumA_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45G
#SBATCH --time=2:00:00

R_SCRIPT="/home/ctools/opt/R-4.4.1/bin/Rscript"

echo "Running BRCA_LumA"
echo "Job ID: ${SLURM_JOB_ID}"

$R_SCRIPT ../R/02_run_tcga_pairedGSEA.R
