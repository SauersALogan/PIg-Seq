#!/bin/bash --login
# =============================================================================
# Script:       14_checkv.sh
# Author:       Kirsten A. L. Y. Lim
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Quality assessment of predicted proviruses using CheckV,
#               followed by filtering to completeness >=50% and exclusion
#               of flagged warnings
# Dependencies: CheckV, CheckV database v1.5
# Input:        Provirus FASTA from geNomad (13_genomad.sh output)
# Output:       Filtered quality summary in CheckV output directory
# Usage:        sbatch 14_checkv.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0
#SBATCH -p multicore
#SBATCH -n 8
#SBATCH -J checkv
#SBATCH -o SLURM_outputs/checkv.out
#SBATCH -e SLURM_errors/checkv.err

# ---- Paths (update for your system) ----
provirus_fasta="/path/to/genomad/output/all_bins_renamed_find_proviruses/all_bins_renamed_provirus.fna"
out_dir="/path/to/checkV/output"
checkv_db="/path/to/checkv-db-v1.5"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate checkV

# ---- Run ----
checkv end_to_end ${provirus_fasta} ${out_dir} \
    -t $SLURM_NTASKS -d ${checkv_db}

## ---- Filter to completeness >=50%, exclude warnings ----
awk -F'\t' 'NR==1 {print; next} \
    $10 ~ /^[0-9.]+$/ && $10 >= 50 && \
    $14 !~ /contig >1.5x longer than expected genome length|high kmer_freq|no viral genes detected/' \
    ${out_dir}/quality_summary.tsv > ${out_dir}/filtered_quality_summary.tsv
