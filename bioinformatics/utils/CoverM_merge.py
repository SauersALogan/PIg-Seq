#!/usr/bin/env python3
"""
Merge CoverM bin abundance TSV files by genome/bin ID.
"""

import pandas as pd
import glob
import argparse
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description='Merge multiple CoverM bin abundance TSV files by Genome ID'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input file pattern with wildcards (e.g., "DNA/CoverM/*_bin_abundance_coverm.tsv")'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output merged TSV file'
    )
    parser.add_argument(
        '--keep-original-headers',
        action='store_true',
        help='Keep original column headers as-is (default: add sample prefix)'
    )

    args = parser.parse_args()

    # Get all matching files
    files = glob.glob(args.input)

    if not files:
        print(f"Error: No files found matching pattern: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files)} files to merge")

    # Read and merge all files
    dfs = []
    for file in files:
        print(f"Reading {file}...")
        df = pd.read_csv(file, sep='\t')

        if not args.keep_original_headers:
            # Get sample name from filename
            sample_name = Path(file).stem.replace('_bin_abundance_coverm', '')
            # Rename columns to include sample name (except 'Genome')
            df.columns = ['Genome'] + [f"{sample_name}_{col}" for col in df.columns[1:]]

        dfs.append(df)

    # Merge all dataframes on 'Genome' column
    print("Merging dataframes...")
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on='Genome', how='outer')

    # Fill NaN values with 0
    merged = merged.fillna(0)

    # Save merged table
    print(f"Writing merged output to {args.output}")
    merged.to_csv(args.output, sep='\t', index=False)
    print(f"Done! Merged {len(merged)} genomes across {len(files)} samples")

if __name__ == '__main__':
    main()
