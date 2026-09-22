#!/bin/bash --login
# =============================================================================
# Script:       05_dereplication.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Dereplication of refined MAGs using dRep2 with fastANI at
#               98% ANI threshold (completeness >90%, contamination <5%)
# Dependencies: dRep2, fastANI
# Input:        Refined bins (.fa) and genome quality CSV
# Output:       Dereplicated bin set in output directory
# Usage:        sbatch 05_dereplication.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p himem
#SBATCH -C icelake
#SBATCH -n 8
#SBATCH -J dereplication
#SBATCH -o SLURM_outputs/dereplication.out
#SBATCH -e SLURM_errors/dereplication.err

# ---- Paths (update for your system) ----
bin_dir="/path/to/refined/bins"
out_dir="/path/to/dereplication/output"
genome_info="/path/to/drep_genome_quality_fa.csv"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate drep2

# ---- Run ----
dRep dereplicate ${out_dir} \
    -g ${bin_dir}/*.fa \
    --S_algorithm fastANI \
    -sa 0.99 \
    -comp 90 \
    -con 5 \
    --genomeInfo ${genome_info}
