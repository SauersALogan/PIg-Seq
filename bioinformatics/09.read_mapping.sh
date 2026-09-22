#!/bin/bash --login
# =============================================================================
# Script:       09_read_mapping.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Map cleaned paired-end reads to dereplicated bin reference
#               using BWA-MEM2
# Dependencies: BWA-MEM2, samtools
# Input:        Cleaned paired-end fastq files per sample
# Output:       SAM alignments in output directory
# Usage:        sbatch 09_read_mapping.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p himem
#SBATCH -C icelake
#SBATCH -n 2
#SBATCH -J read_mapping
#SBATCH -o SLURM_outputs/%A_%a.out
#SBATCH -e SLURM_errors/%A_%a.err
#SBATCH --array=0-23

# ---- Paths (update for your system) ----
read_dir="/path/to/cleaned/reads"
out_dir="/path/to/read/mapping/output"
index="/path/to/bwa-mem2/bin_index"
read_1_path=${read_dir}/${sample_id}_DNA_1.fastq
read_2_path=${read_dir}/${sample_id}_DNA_2.fastq

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate EvoPipe

# ---- Run ----
mapfile -t samples < samples.txt
sample_id=${samples[$SLURM_ARRAY_TASK_ID]}

mkdir -p ${out_dir}

bwa-mem2 mem -t $SLURM_NTASKS ${index} ${read_1_path} ${read_2_path} > ${out_dir}/${sample_id}_mapped.sam
