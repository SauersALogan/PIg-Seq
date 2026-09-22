#!/bin/bash --login
# =============================================================================
# Script:       07_classification.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Taxonomic classification of dereplicated MAGs using GTDB-Tk
# Dependencies: GTDB-Tk, GTDB reference database (release 226)
# Input:        Dereplicated bins (.fa)
# Output:       GTDB-Tk classification results
# Usage:        sbatch 07_classification.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -p himem   # This allows us access to the high memory nodes
#SBATCH -C icelake
#SBATCH -n 4          # (or --ntasks=) Number of cores (2--168 on AMD)
#SBATCH -t 1-0 # This sets the time limit for the job at 24 hours if we did 0-1 it would be one hour
#SBATCH -J Bin_classification
#SBATCH -o SLURM_outputs/Classification.out
#SBATCH -e SLURM_errors/Classification.err

# ---- Paths (update for your system) ----
bin_dir="/path/to/dereplicated/genomes"
out_dir="/path/to/classification/output"
gtdbtk_db="/path/to/GTDB_reference/release226"
mash_db="/path/to/mash/database"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate GTDB

# ---- Run ----
gtdbtk classify_wf --genome_dir ${bin_dir} \
    --out_dir ${out_dir} \
    --extension .fa \
    --mash_db ${mash_db}
