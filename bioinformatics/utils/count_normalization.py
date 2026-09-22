#!/usr/bin/env python3
"""
normalize_bin_counts.py - Normalize and annotate bin-based featureCounts output
"""

import argparse
import pandas as pd
import glob
import os
import sys

def normalize_and_annotate_bins(count_file, gff_file, output_dir):
    """
    Normalize featureCounts output and annotate with product information
    """
    try:
        # Read featureCounts output (skip comment lines)
        counts = pd.read_csv(count_file, delimiter='\t', comment='#')

        # Count column is the last column in featureCounts output
        count_col = counts.columns[-1]
        total_counts = counts[count_col].sum()

        if total_counts == 0:
            print(f"Warning: No counts found in {count_file}")
            return None

        # Calculate normalizations
        counts['CPM'] = (counts[count_col] / total_counts) * 1_000_000
        counts['RPK'] = counts[count_col] / (counts['Length'] / 1000)
        total_RPK = counts['RPK'].sum()

        if total_RPK > 0:
            counts['TPM'] = (counts['RPK'] / total_RPK) * 1_000_000
        else:
            counts['TPM'] = 0

        # Read GFF and extract annotations
        gff = pd.read_csv(gff_file, delimiter='\t', comment='#', header=None)
        cds_df = gff[gff[2] == 'CDS'].copy()

        # For bin-based system, extract ID instead of locus_tag
        cds_df['gene_id'] = cds_df[8].str.extract(r'ID=([^;]+)')
        cds_df['product'] = cds_df[8].str.extract(r'product=([^;]+)')
        cds_df['product'] = cds_df['product'].fillna('hypothetical protein')

        # Create annotation mapping
        annotations = dict(zip(cds_df['gene_id'], cds_df['product']))
        counts['product'] = counts['Geneid'].map(annotations)
        counts['product'] = counts['product'].fillna('hypothetical protein')

        # Create output filename
        count_base = os.path.basename(count_file)
        count_name = os.path.splitext(count_base)[0]
        output_name = count_name + "_normalized.txt"
        output_path = os.path.join(output_dir, output_name)

        # Save output
        counts.to_csv(output_path, sep='\t', index=False)
        print(f"Normalized counts saved to: {output_path}")
        return output_path

    except Exception as e:
        print(f"Error processing {count_file}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Normalize and annotate bin-based featureCounts output")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--count_files", nargs="+", required=True, help="Count files or glob patterns")
    parser.add_argument("--gff_file", required=True, help="Single merged GFF file for annotation")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Expand glob patterns for count files
    count_files = []
    for pattern in args.count_files:
        expanded = glob.glob(pattern)
        if expanded:
            count_files.extend(expanded)
        else:
            # If no glob match, treat as literal filename
            if os.path.exists(pattern):
                count_files.append(pattern)
            else:
                print(f"Warning: No files found matching pattern: {pattern}")

    count_files.sort()

    if not count_files:
        print("Error: No count files found")
        sys.exit(1)

    if not os.path.exists(args.gff_file):
        print(f"Error: GFF file not found: {args.gff_file}")
        sys.exit(1)

    if args.verbose:
        print(f"Found {len(count_files)} count files")
        print(f"Using GFF file: {args.gff_file}")
        print(f"Output directory: {args.output}")

    # Process each count file
    successful = 0
    failed = 0

    for count_file in count_files:
        if args.verbose:
            print(f"Processing: {count_file}")

        result = normalize_and_annotate_bins(count_file, args.gff_file, args.output)

        if result:
            successful += 1
        else:
            failed += 1

    print(f"\nProcessing complete:")
    print(f"Successfully processed: {successful} files")
    print(f"Failed: {failed} files")

if __name__ == "__main__":
    main()
