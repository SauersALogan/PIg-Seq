#!/bin/bash --login
# =============================================================================
# Script:       12_subsample_reads.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Subsample cleaned reads to equal depth across all samples
#               using seqtk, then map to dereplicated bin reference
# Dependencies: seqtk, bwa-mem2, samtools
# Input:        Cleaned paired-end fastq files per sample
# Output:       Subsampled BAMs in subsampled_bam_files/
# Usage:        sbatch 12_subsample_reads.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 1-0
#SBATCH -p multicore
#SBATCH -n 4
#SBATCH -J subsample_map
#SBATCH -o SLURM_outputs/%A_%a.out
#SBATCH -e SLURM_errors/%A_%a.err
#SBATCH --array=0-23

# ---- Paths (update for your system) ----
read_dir="/path/to/cleaned/reads"
out_dir="/path/to/subsampled/output"
index="/path/to/bwa-mem2/index"
seed=42
target=3015000

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate EvoPipe

# ---- Run ----
mapfile -t samples < samples.txt
sample_id=${samples[$SLURM_ARRAY_TASK_ID]}

mkdir -p ${out_dir}

seqtk sample -s ${seed} ${read_dir}/${sample_id}_DNA_1.fastq ${target} > ${out_dir}/${sample_id}_sub_1.fastq
seqtk sample -s ${seed} ${read_dir}/${sample_id}_DNA_2.fastq ${target} > ${out_dir}/${sample_id}_sub_2.fastq

bwa-mem2 mem -t $SLURM_NTASKS ${index} \
    ${out_dir}/${sample_id}_sub_1.fastq \
    ${out_dir}/${sample_id}_sub_2.fastq | \
    samtools view -bS -q 20 | \
    samtools sort -o ${out_dir}/${sample_id}_mapped.bam

samtools index ${out_dir}/${sample_id}_mapped.bam
