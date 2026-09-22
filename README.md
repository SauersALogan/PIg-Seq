# PIg-Seq
Scripts for "Paired sequencing of IgA-bound bacteria reveals widespread associations between adaptive immunity and gut microbiome gene expression."

Processed data available at https://doi.org/10.48420/33951214. Raw reads deposited at ENA under PRJEB127118.

# Repository structure
Relevant bioinformatics scripts and utilities can be found in the bioinformatics/ folder. Scripts for statistical analysis and related functions can be found in the statistics/ folder.

# Bioinformatics pipeline
Bioinformatics scripts 01–08 must be run before the phage scripts.

Script 09 (read_mapping.sh) is used for mapping all reads and subsampled reads to the dereplicated bins, along with mapping reads to phage genomes.

Scripts 10 and 11 were used for all reads, subsampled reads, and phage-specific reads.

Script 12 subsamples cleaned reads to equal depth across samples using seqtk before mapping (script 09) and counting (scripts 10–11) on the subsampled set.

Scripts 13–18 cover the phage analysis pipeline: provirus detection (geNomad), quality filtering (CheckV), viral dereplication at 95% ANI / 85% target coverage, and annotation (Pharokka, phold, Phynteny).

The utility scripts in bioinformatics/utils/ accept featureCounts and CoverM output files and normalize and merge sample results into unified .tsv count matrices.

# Environment
Environment setup in shell scripts is specific to the University of Manchester CSF3 cluster. Users on other systems will need to install the listed dependencies independently.

# Statistical analysis

Shared R functions sourced by all analysis scripts are in statistics/functions/ 

Analysis scripts 01–16 in statistics/ are numbered loosely by manuscript figure order. Each script declares its dependencies, input files, and figure output in its header. Script 07a is a Python script producing intermediate CSVs consumed by 07b.

The data/ and figures/ directories are gitignored. Download processed data from the Figshare DOI above and place files in data/.
