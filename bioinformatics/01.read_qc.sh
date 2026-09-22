#!/bin/bash --login
# =============================================================================
# Script:       01_read_qc.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Adapter trimming (Trimmomatic) and host read removal (bowtie2)
#               for paired-end metagenomic and metatranscriptomic reads
# Dependencies: Trimmomatic 0.39, bowtie2, metawrap
# Input:        Raw paired-end .fq files per sample
# Output:       Host-filtered paired reads in DNA/readqc/<sample_id>/
# Usage:        sbatch 01_read_qc.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p multicore
#SBATCH -n 4
#SBATCH -J read_qc-pipeline
#SBATCH -o SLURM_outputs/%A_%a.out
#SBATCH -e SLURM_errors/%A_%a.err
#SBATCH --array=0-23

# ---- Database paths (cluster-specific, update for your system) ----
genome_index="/path/to/mouse/bowtie2/index"
read_dir="/path/to/raw/reads"
outdir="/path/to/output"

# ---- Environment ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
module load apps/binapps/trimmomatic/0.39
mamba activate metawrap

# ---- Run ----
mapfile -t samples < samples.txt
sample_id=${samples[$SLURM_ARRAY_TASK_ID]}

bash bioinformatics/utils/qc_pipeline.sh \
    --outdir=${outdir}/${sample_id} \
    --DNA1=${read_dir}/${sample_id}_DNA_1.fq \
    --DNA2=${read_dir}/${sample_id}_DNA_2.fq \
    --DNAOnly --threads=$SLURM_NTASKS \
    --genome_index=${genome_index}
