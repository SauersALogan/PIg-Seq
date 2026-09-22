#!/bin/bash --login
# =============================================================================
# Script:       06_checkm2.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Quality assessment of refined MAGs using CheckM2
# Dependencies: CheckM2, UniRef100 KO diamond database
# Input:        Refined bins (.fa)
# Output:       CheckM2 quality report in output directory
# Usage:        sbatch 06_checkm2.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p himem
#SBATCH -C cascadelake
#SBATCH -n 6
#SBATCH -J CheckM2
#SBATCH -o SLURM_outputs/CheckM2.out
#SBATCH -e SLURM_errors/CheckM2.err

# ---- Paths (update for your system) ----
bin_dir="/path/to/refined/bins"
out_dir="/path/to/checkm2/results"
checkm2_db="/path/to/CheckM2_database/uniref100.KO.1.dmnd"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate UPNGS-refinement

# ---- Run ----
checkm2 predict --threads $SLURM_NTASKS \
    --input ${bin_dir} \
    --output-directory ${out_dir} \
    --extension fa \
    --database_path ${checkm2_db} \
    --force
