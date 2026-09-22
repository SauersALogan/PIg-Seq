#!/bin/bash --login
# =============================================================================
# Script:       02_assembly.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Metagenomic assembly using metaSPAdes via metaWRAP
# Dependencies: metaWRAP, SPAdes
# Input:        QC-filtered paired-end .fastq files per sample
# Output:       Assembled contigs in DNA/assembly/raw_v2_<sample_id>/
# Usage:        sbatch 02_assembly.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p himem
#SBATCH -C icelake
#SBATCH -n 12
#SBATCH -J assembly_metaspades_only
#SBATCH -o SLURM_outputs/%A_%a.out
#SBATCH -e SLURM_errors/%A_%a.err
#SBATCH --array=0-23

# ---- Paths (update for your system) ----
read_dir="/path/to/qc/filtered/reads"
outdir="/path/to/assembly/output"
read_1_path=${read_dir}/${sample_id}_DNA_1.fastq
read_2_path=${read_dir}/${sample_id}_DNA_2.fastq

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate metawrap

# ---- Run ----
mapfile -t samples < samples.txt
sample_id=${samples[$SLURM_ARRAY_TASK_ID]}

metawrap assembly -t $SLURM_NTASKS -o ${outdir}/${sample_id} \
    -1 ${read_1_path} \
    -2 ${read_2_path} \
    -m 800 --metaspades
