#!/bin/bash --login
# =============================================================================
# Script:       16_pharokka.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Base annotation of dereplicated provirus genomes using Pharokka
# Dependencies: Pharokka, Pharokka database
# Input:        Filtered dereplicated provirus FASTA (15_viral_dereplication.sh
#               output)
# Output:       Pharokka annotations (GBK, GFF, FAA, FNA) in output directory
# Usage:        sbatch 16_pharokka.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p multicore_small
#SBATCH -n 4
#SBATCH -J pharokka
#SBATCH -o SLURM_outputs/pharokka
#SBATCH -e SLURM_errors/pharokka

# ---- Paths (update for your system) ----
provirus_fasta="/path/to/filtered_provirus.fna"
out_dir="/path/to/pharokka/output"
pharokka_db="/path/to/pharokka/db"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate pharokka_env

# ---- Run ----
pharokka.py -i ${provirus_fasta} \
    -o ${out_dir} \
    -d ${pharokka_db} \
    -t $SLURM_NTASKS \
    --force
