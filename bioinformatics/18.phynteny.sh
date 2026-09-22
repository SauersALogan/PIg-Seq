#!/bin/bash --login
# =============================================================================
# Script:       18_phynteny.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Synteny-based functional prediction of provirus genes using
#               Phynteny transformer, layered on phold output
# Dependencies: Phynteny, Phynteny model
# Input:        phold GBK output (17_phold.sh output)
# Output:       Phynteny annotations in output directory
# Usage:        sbatch 18_phynteny.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p gpuL
#SBATCH -G 2
#SBATCH -n 2
#SBATCH -J phynteny
#SBATCH -o SLURM_outputs/phynteny.out
#SBATCH -e SLURM_errors/phynteny.err

# ---- Paths (update for your system) ----
phold_gbk="/path/to/phold/output/phage_annotate.gbk"
out_dir="/path/to/phynteny/output"
phynteny_model="/path/to/phynteny/model"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate phynteny

# ---- Run ----
phynteny_transformer ${phold_gbk} \
    -o ${out_dir} \
    -m ${phynteny_model}
