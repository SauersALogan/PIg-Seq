#!/bin/bash --login
# =============================================================================
# Script:       11_featureCounts.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Gene-level read counting using featureCounts on mapped reads
#               against Bakta GFF annotations. Counts only properly paired,
#               concordant fragments with MAPQ ≥20, excluding chimeric reads.
# Dependencies: Subread (featureCounts)
# Input:        SAM alignments and filtered GFF annotation
# Output:       Per-sample gene count tables
# Usage:        sbatch 11_featureCounts.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p multicore_small
#SBATCH -n 2
#SBATCH -J featureCounts
#SBATCH -o SLURM_outputs/%A_%a.out
#SBATCH -e SLURM_errors/%A_%a.err
#SBATCH --array=0-23

# ---- Paths (update for your system) ----
sam_dir="/path/to/read/mapping/output"
out_dir="/path/to/featureCounts/output"
gff="/path/to/filtered.gff"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate EvoPipe

# ---- Run ----
mapfile -t samples < samples.txt
sample_id=${samples[$SLURM_ARRAY_TASK_ID]}

featureCounts -T 1 -t CDS -g ID -a ${gff} \
    -o ${out_dir}/${sample_id}_counts.txt -p -B -C -Q 20 \
    ${sam_dir}/${sample_id}_mapped.sam
