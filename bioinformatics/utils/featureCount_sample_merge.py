#!/usr/bin/env python3
"""
merge_normalized_counts.py - Merge normalized count files into sample matrices
"""

import argparse
import pandas as pd
import glob
import os
import sys

def extract_sample_name(filename):
    """Extract sample name from filename"""
    basename = os.path.basename(filename)
    # Remove _counts_normalized.txt suffix
    sample_name = basename.replace('_counts_normalized.txt', '')
    return sample_name

def merge_normalized_files(input_files, output_file, metric='TPM'):
    """
    Merge normalized count files into a matrix where samples are columns

    Parameters:
    - input_files: List of normalized count files
    - output_file: Output file path
    - metric: Which metric to extract (CPM, RPK, TPM, or raw counts)
    """

    if not input_files:
        print("Error: No input files provided")
        return None

    all_data = []
    sample_names = []

    # Read each file and extract the desired metric
    for file_path in input_files:
        try:
            df = pd.read_csv(file_path, sep='\t')
            sample_name = extract_sample_name(file_path)
            sample_names.append(sample_name)

            # Determine which column to use
            if metric in ['CPM', 'RPK', 'TPM']:
                if metric not in df.columns:
                    print(f"Warning: {metric} column not found in {file_path}")
                    continue
                count_col = metric
            else:  # assume raw counts (last .sam column)
                sam_cols = [col for col in df.columns if '.sam' in col]
                if not sam_cols:
                    print(f"Warning: No .sam column found in {file_path}")
                    continue
                count_col = sam_cols[0]

            # Create sample dataframe with gene info and counts
            sample_df = df[['Geneid', 'Chr', 'Length', count_col, 'product']].copy()
            sample_df = sample_df.rename(columns={count_col: sample_name})

            all_data.append(sample_df)
            print(f"Processed {sample_name}: {len(sample_df)} genes")

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue

    if not all_data:
        print("Error: No files successfully processed")
        return None

    # Merge all samples on gene information
    merged_df = all_data[0].copy()

    for i, sample_df in enumerate(all_data[1:], 1):
        # Merge on gene identifiers, keeping Chr, Length, product from first file
        merged_df = merged_df.merge(
            sample_df[['Geneid', sample_names[i]]],
            on='Geneid',
            how='outer'
        )

    # Fill missing values with 0
    count_columns = sample_names
    merged_df[count_columns] = merged_df[count_columns].fillna(0)

    # Reorder columns: Geneid, Chr, Length, product, then all sample columns
    column_order = ['Geneid', 'Chr', 'Length', 'product'] + count_columns
    merged_df = merged_df[column_order]

    # Save merged matrix
    merged_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nMerged matrix saved to: {output_file}")
    print(f"Matrix dimensions: {len(merged_df)} genes x {len(count_columns)} samples")

    return output_file

def main():
    parser = argparse.ArgumentParser(description="Merge normalized count files into sample matrix")
    parser.add_argument("--input", nargs="+", required=True,
                       help="Normalized count files or glob patterns")
    parser.add_argument("--output", required=True, help="Output file path (.txt)")
    parser.add_argument("--metric", choices=['CPM', 'RPK', 'TPM', 'raw'], default='TPM',
                       help="Which metric to extract (default: TPM)")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")

    args = parser.parse_args()

    # Expand glob patterns
    input_files = []
    for pattern in args.input:
        expanded = glob.glob(pattern)
        if expanded:
            input_files.extend(expanded)
        else:
            if os.path.exists(pattern):
                input_files.append(pattern)
            else:
                print(f"Warning: No files found matching pattern: {pattern}")

    input_files.sort()

    if not input_files:
        print("Error: No input files found")
        sys.exit(1)

    if args.verbose:
        print(f"Found {len(input_files)} input files")
        print(f"Extracting {args.metric} values")
        print(f"Output file: {args.output}")

    # Merge files
    result = merge_normalized_files(input_files, args.output, args.metric)

    if result:
        print("\nMerging completed successfully!")
    else:
        print("Error: Merging failed")
        sys.exit(1)

if __name__ == "__main__":
    main()
