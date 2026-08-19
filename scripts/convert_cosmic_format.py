#!/usr/bin/env python3

from argparse import ArgumentParser
import pandas as pd

def parse_arguments():
    parser = ArgumentParser(description="Convert COSMIC signature file to different formats")
    parser.add_argument("-i", "--input", help="Input COSMIC signature file", required=True)
    parser.add_argument("-o", "--output", help="Output file path", required=True)
    parser.add_argument("-f", "--format", help="Output format (csv_no_index, csv_quoted, csv_with_index, tab_with_index)", required=True)
    parser.add_argument("-r", "--relevant-signatures", help="Comma-separated list of relevant signatures to filter TO (include only these)", default=None)
    parser.add_argument("-e", "--exclude-signatures", help="Comma-separated list of signatures to EXCLUDE (out-of-reference)", default=None)
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Read base COSMIC file
    cosmic_df = pd.read_csv(args.input, sep='\t', index_col='Type')
    
    # Filter to relevant signatures if specified (INCLUDE only these)
    if args.relevant_signatures:
        relevant_sigs = [sig.strip() for sig in args.relevant_signatures.split(',')]
        available_sigs = [sig for sig in relevant_sigs if sig in cosmic_df.columns]
        cosmic_df = cosmic_df[available_sigs]
    
    # Exclude out-of-reference signatures if specified (REMOVE these)
    if args.exclude_signatures:
        exclude_sigs = [sig.strip() for sig in args.exclude_signatures.split(',')]
        remaining_sigs = [sig for sig in cosmic_df.columns if sig not in exclude_sigs]
        cosmic_df = cosmic_df[remaining_sigs]
    
    # Apply format-specific settings
    if args.format == "csv_no_index":
        # Comma-separated without row names (MutationalPatterns, mmsig, YAPSA)
        cosmic_df.to_csv(args.output, sep=',', index=False, header=True)
    elif args.format == "csv_quoted":
        # Comma-separated with quotes around all values (sigLASSO)
        cosmic_df.to_csv(args.output, sep=',', index=True, header=True, quoting=1)
    elif args.format == "csv_with_index":
        # Comma-separated with row names (sigfit)
        cosmic_df.to_csv(args.output, sep=',', index=True, header=True)
    elif args.format == "tab_with_index":
        # Tab-separated with row names (default, deconstructSigs, sigminer, etc.)
        cosmic_df.to_csv(args.output, sep='\t', index=True, header=True)
    else:
        raise ValueError(f"Unknown format: {args.format}")

if __name__ == "__main__":
    main()
