#!/bin/bash --login
# =============================================================================
# Script:       04_bin_refinement.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Bin refinement using metaWRAP bin_refinement, consolidating
#               MetaBAT2, MaxBin2, and CONCOCT bin sets (completeness >90%,
#               contamination <5%)
# Dependencies: metaWRAP
# Input:        Binned MAGs from MetaBAT2, MaxBin2, and CONCOCT per sample
# Output:       Refined bins in DNA/refined_bins/<sample_id>/
# Usage:        sbatch 04_bin_refinement.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p himem
#SBATCH -C cascadelake
#SBATCH -n 2
#SBATCH -J refinement
#SBATCH -o SLURM_outputs/%A_%a.out
#SBATCH -e SLURM_errors/%A_%a.err
#SBATCH --array=0-23

# ---- Paths (update for your system) ----
bin_dir="/path/to/binning/output"
out_dir="/path/to/refined/bins"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate metawrap

# ---- Run ----
mapfile -t samples < samples.txt
sample_id=${samples[$SLURM_ARRAY_TASK_ID]}

metawrap bin_refinement -o ${out_dir}/${sample_id} -t $SLURM_NTASKS \
    -A ${bin_dir}/${sample_id}/maxbin2_bins \
    -B ${bin_dir}/${sample_id}/metabat2_bins \
    -C ${bin_dir}/${sample_id}/concoct_bins \
    -c 90 -x 5
