#!/bin/bash --login
# =============================================================================
# Script:       17_phold.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Protein structure-based annotation of provirus genomes using
#               phold with GPU-accelerated foldseek, layered on Pharokka output
# Dependencies: phold, phold database, GPU node
# Input:        Pharokka GBK output (16_pharokka.sh output)
# Output:       phold annotations in output directory
# Usage:        sbatch 17_phold.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p gpuL
#SBATCH -G 2
#SBATCH -n 2
#SBATCH -J phold
#SBATCH -o SLURM_outputs/phold
#SBATCH -e SLURM_errors/phold

# ---- Paths (update for your system) ----
pharokka_gbk="/path/to/pharokka/output/pharokka.gbk"
out_dir="/path/to/phold/output"
phold_db="/path/to/phold/db"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate phold_new

# ---- Run ----
phold run -i ${pharokka_gbk} \
    -o ${out_dir} \
    -t $SLURM_NTASKS \
    --prefix phage_annotate \
    -d ${phold_db} \
    --foldseek_gpu \
    --force
