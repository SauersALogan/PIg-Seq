#!/bin/bash --login
# =============================================================================
# Script:       08_annotation.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Gene annotation of dereplicated MAGs using Bakta, with
#               bin-specific locus tag prefixes
# Dependencies: Bakta, Bakta database
# Input:        Dereplicated bins (.fa), bins.txt sample list
# Output:       Annotated genomes (GFF, FAA, FNA) per bin
# Usage:        sbatch 08_annotation.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p multicore
#SBATCH -n 6
#SBATCH -J bin_annotation
#SBATCH -o SLURM_outputs/%A_%a.out
#SBATCH -e SLURM_errors/%A_%a.err
#SBATCH --array=0-89

# ---- Paths (update for your system) ----
bin_dir="/path/to/dereplicated/genomes"
out_dir="/path/to/bin/annotations"
bakta_db="/path/to/bakta/db"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate bakta

# ---- Run ----
mapfile -t samples < bins.txt
sample_id=${samples[$SLURM_ARRAY_TASK_ID]}

bakta --db ${bakta_db} \
    --output ${out_dir}/${sample_id} \
    --locus "Bin.${SLURM_ARRAY_TASK_ID}" \
    --locus-tag "Bin.${SLURM_ARRAY_TASK_ID}.locus" \
    --prefix "${sample_id}" \
    --threads $SLURM_NTASKS \
    --verbose \
    ${bin_dir}/${sample_id}.fa
