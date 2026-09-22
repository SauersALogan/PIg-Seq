#!/bin/bash --login
# =============================================================================
# Script:       10_coverm_abundance.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Convert SAM to filtered sorted BAM (MAPQ ≥20), then estimate
#               bin-level TPM abundances using CoverM
# Dependencies: samtools, CoverM
# Input:        SAM alignments and genome definitions file
# Output:       Sorted BAMs and per-sample CoverM TPM abundance tables
# Usage:        sbatch 10_coverm_abundance.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p multicore_small
#SBATCH -n 2
#SBATCH -J CoverM
#SBATCH -o SLURM_outputs/CoverM_%A_%a.out
#SBATCH -e SLURM_errors/CoverM_%A_%a.err
#SBATCH --array=0-23

# ---- Paths (update for your system) ----
sam_dir="/path/to/read/mapping/output"
bam_dir="/path/to/bam/output"
out_dir="/path/to/coverm/abundance/output"
genome_def="/path/to/genome_definitions.tsv"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate CoverM

# ---- Run ----
mapfile -t samples < samples.txt
sample_id=${samples[$SLURM_ARRAY_TASK_ID]}

mkdir -p ${bam_dir} ${out_dir}

samtools view -q 20 -bS ${sam_dir}/${sample_id}_mapped.sam \
    | samtools sort -@ ${SLURM_NTASKS} -o ${bam_dir}/${sample_id}.bam

samtools index ${bam_dir}/${sample_id}.bam

coverm genome \
    --bam-files ${bam_dir}/${sample_id}.bam \
    --genome-definition ${genome_def} \
    -t $SLURM_NTASKS -m tpm -o ${out_dir}/${sample_id}_bin_abundance_coverm.tsv
