#!/bin/bash --login
# =============================================================================
# Script:       13_genomad.sh
# Author:       Kirsten A. L. Y. Lim
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Identify proviruses in dereplicated MAGs using geNomad
# Dependencies: geNomad, geNomad database
# Input:        Dereplicated bin FASTA (all_bins_renamed.fa)
# Output:       geNomad provirus predictions in output directory
# Usage:        sbatch 13_genomad.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0
#SBATCH -p multicore
#SBATCH -n 4
#SBATCH -J genomad
#SBATCH -o SLURM_outputs/genomad.out
#SBATCH -e SLURM_errors/genomad.err

# ---- Paths (update for your system) ----
bin_fasta="/path/to/dereplicated/all_bins_renamed.fa"
out_dir="/path/to/genomad/output"
genomad_db="/path/to/genomad_db"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate genomad

# ---- Run ----
genomad end-to-end --splits 10 ${bin_fasta} ${out_dir} ${genomad_db}
