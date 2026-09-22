#!/bin/bash --login
# =============================================================================
# Script:       15_viral_dereplication.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Dereplication of predicted proviruses at 95% ANI / 85% target
#               coverage following MIUViG species-level thresholds
#               (Roux et al. 2019, Nature Biotechnology). Uses CheckV's
#               anicalc/aniclust approach via all-vs-all BLASTn.
# Dependencies: seqkit, BLASTn, CheckV scripts: anicalc.py, aniclust.py
# Input:        Provirus FASTA from geNomad, quality CSV from CheckV
#               (14_checkv.sh output)
# Output:       Dereplicated viral clusters (viral_clusters.tsv)
# Usage:        sbatch 15_viral_dereplication.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p multicore_small
#SBATCH -n 4
#SBATCH -J dereplication
#SBATCH -o SLURM_outputs/viral_dereplication.out
#SBATCH -e SLURM_errors/viral_dereplication.err

# ---- Paths (update for your system) ----
provirus_fasta="/path/to/all_bins_renamed_provirus.fna"
quality_csv="/path/to/viral_genome_quality.csv"
out_dir="/path/to/viral_analysis/output"
checkv_scripts="/path/to/checkv/scripts"

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate EvoPipe

## ---- Filter to completeness >50% ----
awk -F',' 'NR>1 && $2>50 {gsub(/\.fna$/, "", $1); print $1}' \
    ${quality_csv} > ${out_dir}/selected_genomes.txt

wc -l ${out_dir}/selected_genomes.txt

seqkit grep -f ${out_dir}/selected_genomes.txt \
    ${provirus_fasta} > ${out_dir}/selected_provirus.fna

## ---- All-vs-all BLASTn ----
blastn -query ${out_dir}/selected_provirus.fna \
    -subject ${out_dir}/selected_provirus.fna \
    -outfmt "6 std qlen slen" \
    -out ${out_dir}/viral_blastn.tsv \
    -perc_identity 90 \
    -max_target_seqs 10000

## ---- ANI calculation and clustering (95% ANI, 85% target coverage) ----
python ${checkv_scripts}/anicalc.py -i ${out_dir}/viral_blastn.tsv -o ${out_dir}/viral_ani.tsv

python ${checkv_scripts}/aniclust.py --fna ${out_dir}/selected_provirus.fna \
    --ani ${out_dir}/viral_ani.tsv \
    --out ${out_dir}/viral_clusters.tsv \
    --min_ani 95 \
    --min_qcov 0 \
    --min_tcov 85
