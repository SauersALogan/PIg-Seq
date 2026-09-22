#!/bin/bash --login
# =============================================================================
# Script:       03_binning.sh
# Author:       Logan A. Sauers
# Repository:   github.com/SauersALogan/PIg-Seq
# Description:  Metagenomic binning using MetaBAT2 and MaxBin2 (via metaWRAP)
#               followed by CONCOCT binning in a separate environment
# Dependencies: metaWRAP (MetaBAT2, MaxBin2), CONCOCT, samtools
# Input:        Assembled contigs and QC-filtered paired reads per sample
# Output:       Binned MAGs in DNA/binning/<sample_id>/
# Usage:        sbatch 03_binning.sh
# =============================================================================

# ---- SLURM configuration ----
#SBATCH -t 4-0           # days-hours
#SBATCH -p multicore
#SBATCH -n 4
#SBATCH -J binning
#SBATCH -o SLURM_outputs/%A_%a.out
#SBATCH -e SLURM_errors/%A_%a.err
#SBATCH --array=0-23

# ---- Paths (update for your system) ----
assembly_dir="/path/to/assemblies"
read_dir="/path/to/cleaned/reads"
bin_dir="/path/to/binning/output"
assembly=${assembly}/${sample_id}/final_assembly.fasta
outdir=${bin_dir}/${sample_id}

# ---- Environment (CSF3-specific) ----
module use -a /mnt/bmh01-rds/coyte_lab_sequencing/software/common/modulesfiles
module load apps/binapps/anaconda3/2023.03
module load conda_envs
mamba activate metawrap

# ---- Run ----
mapfile -t samples < samples.txt
sample_id=${samples[$SLURM_ARRAY_TASK_ID]}

metawrap binning -o ${outdir} \
        -a ${assembly} \
        -t $SLURM_NTASKS --metabat2 --maxbin2 \
	${read_dir}/${sample_id}_DNA_*

# ---- Switch to CONCOCT environment ----
mamba deactivate
mamba activate UPNGS-binning


outdir="DNA/binning/${sample_id}"
assembly="DNA/assembly/${sample_id}/final_assembly.fasta"

samtools index $outdir/work_files/${sample_id}_DNA.bam

# ---- CONCOCT binning ----
if [ ! -d "$outdir/concoct_data" ]; then
	echo ""
	echo "The concoct folder does not exist, creating it now";
	mkdir -p $outdir/concoct_data;
	echo "";
	echo "Run concoct binning"
		cut_up_fasta.py $assembly -c 10000 -o 0 --merge_last \
			-b $outdir/concoct_contigs.bed > $outdir/concoct_10K.fa
		concoct_coverage_table.py $outdir/concoct_contigs.bed \
			$outdir/work_files/${sample_id}_DNA.bam > $outdir/concoct_depth.txt
		concoct --composition_file $outdir/concoct_10K.fa \
			--coverage_file $outdir/concoct_depth.txt -b $outdir/concoct_data/ \
			--threads $SLURM_NTASKS
	echo "Remerge the contigs and split into bins"

elif [ ! "$(ls -A "$outdir/concoct_data")" ]; then
	echo "Concoct folder exists but is empty, I will run concoct now"
	cut_up_fasta.py $assembly -c 10000 -o 0 --merge_last \
			-b $outdir/concoct_contigs.bed > $outdir/concoct_10K.fa
		concoct_coverage_table.py $outdir/concoct_contigs.bed \
			$outdir/work_files/${sample_id}_DNA.bam > $outdir/concoct_depth.txt
		concoct --composition_file $outdir/concoct_10K.fa \
			--coverage_file $outdir/concoct_depth.txt -b $outdir/concoct_data/ \
			--threads $SLURM_NTASKS
else
	echo ""
	echo "Concoct folder exists and is not empty, I dislike overwriting directories"
	echo "I will skip concoct and proceed"
	echo ""
fi

if 	[ ! -d "$outdir/concoct_bins" ]; then
	echo "Concoct bin folder does not exist, creating and splitting bins"
	mkdir $outdir/concoct_bins
	merge_cutup_clustering.py $outdir/concoct_data/clustering_gt1000.csv > $outdir/clustering_merged.csv
	extract_fasta_bins.py $assembly $outdir/clustering_merged.csv \
	--output_path $outdir/concoct_bins
elif [ ! "$(ls -A "$outdir/concoct_bins")" ]; then
	echo "Concoct bin folder exists, but is empty. Proceeding with bin splitting"
	echo "Remerge the contigs and split into bins"
	merge_cutup_clustering.py $outdir/concoct_data/clustering_gt1000.csv > $outdir/clustering_merged.csv
	extract_fasta_bins.py $assembly $outdir/clustering_merged.csv \
	--output_path $outdir/concoct_bins
else
	echo "Concoct bin folder already exists, I dislike overwritting will skip"
	echo "Bin splitting, if you have reached this message in error check logs"
fi
