#!/usr/bin/env python3
"""
Subset a real PCAWG mutational catalog to one cancer type and write it out in the
formats needed for tool fitting. No downsampling: every sample keeps its native
mutation count. Samples are selected via the "CancerType::SampleID" column-naming
convention used in the PCAWG input files.
"""

from argparse import ArgumentParser
import pandas as pd
import os


def parse_arguments():
    parser = ArgumentParser(description="Subset a real catalog to one cancer type, no downsampling")
    parser.add_argument("--input_catalog", required=True,
                        help="Real mutational catalog (contexts x samples, columns named 'CancerType::SampleID')")
    parser.add_argument("--cancer_type", required=True,
                        help="Cancer type to select (must match the prefix before '::' in sample columns)")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--output_formats", nargs='+',
                        default=["deconstructSigs", "MutationalPatterns"])
    parser.add_argument("--num_muts_label", default="none",
                        help="Path label only — no downsampling is performed")
    return parser.parse_args()


def select_cancer_type(real_data, cancer_type):
    """Return only the columns whose 'CancerType::SampleID' prefix matches cancer_type."""
    prefix = f"{cancer_type}::"
    selected = [c for c in real_data.columns if c.startswith(prefix)]
    if not selected:
        available = sorted({c.split("::")[0] for c in real_data.columns if "::" in c})
        raise ValueError(f"No samples found for cancer_type={cancer_type!r}. Available: {available}")
    return real_data[selected]


def save_in_format(counts, output_path, format_name):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    if format_name == "sigfit":
        counts.to_csv(output_path, sep=',', index=True, header=True)
    elif format_name == "sigLASSO":
        counts.to_csv(output_path, sep=',', index=True, header=True, quoting=1)
    else:
        counts.to_csv(output_path, sep='\t', index=True, header=True)


def main():
    args = parse_arguments()

    real_data = pd.read_csv(args.input_catalog, sep='\t', index_col='Type', comment='#')
    counts = select_cancer_type(real_data, args.cancer_type)
    print(f"Using {counts.shape[1]} real {args.cancer_type} samples (native mutation counts, no downsampling)")

    stem = f"{args.cancer_type}-w_0-{args.num_muts_label}"
    for format_name in args.output_formats:
        output_path = os.path.join(args.output_dir, f"counts_{stem}_{format_name}.dat")
        save_in_format(counts, output_path, format_name)
        print(f"Created {output_path}")


if __name__ == "__main__":
    main()
