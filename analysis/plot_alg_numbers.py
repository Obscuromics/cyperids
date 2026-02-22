#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description="Plot number of markers (rows) per category in the 3rd column."
    )
    parser.add_argument("tsv_file", help="Input TSV file")
    parser.add_argument("-o", "--out", help="Output image file (optional)",
                        default=None)
    args = parser.parse_args()

    # Load TSV
    df = pd.read_csv(args.tsv_file, sep="\t")

    # Select the third column regardless of name
    third_col = df.columns[2]

    # Count entries
    counts = df[third_col].value_counts().sort_index()

    # If no ALGs, skip plotting cleanly
    if counts.empty:
        print(f"WARNING: no ALGs found in {args.tsv_file}, skipping plot")
        return

    # Plot
    plt.figure(figsize=(6,4))
    counts.plot(kind="bar")
    plt.xlabel(f"ALG ({third_col})")
    plt.ylabel("Number of markers")
    plt.title(f"Marker count per {third_col}")
    plt.tight_layout()

    # Save or show
    if args.out:
        plt.savefig(args.out, dpi=300)
    else:
        plt.show()

if __name__ == "__main__":
    main()
