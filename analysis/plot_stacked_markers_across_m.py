#!/usr/bin/env python3

import argparse
import glob
import os
import re
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt

from itertools import cycle

def main():
    parser = argparse.ArgumentParser(
        description="Stacked bar plot of marker counts per ALG across m"
    )
    parser.add_argument("--prefix", required=True,
                        help="Prefix, e.g. cyperaceae_topology2")
    parser.add_argument("--node", required=True,
                        help="Node ID, e.g. n2")
    parser.add_argument("--out", default=None,
                        help="Output prefix (optional)")

    args = parser.parse_args()

    prefix = args.prefix
    node = args.node
    out = args.out or f"{node}_{prefix}_marker_stacked_by_m"

    OKABE_ITO = [
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#000000",
    ]

    rows = []
    dir_re = re.compile(r"m(\d+)_output")

    # --------------------------------------------------
    # Read all m*_output directories
    # --------------------------------------------------

    for d in sorted(glob.glob("m*_output")):
        m_match = dir_re.match(d)
        if not m_match:
            continue

        m = int(m_match.group(1))
        tsv = os.path.join(
            d,
            f"m{m}.buscona_{node}.{prefix}.tsv"
        )

        if not os.path.exists(tsv):
            print(f"WARNING: missing {tsv}, skipping")
            continue

        df = pd.read_csv(tsv, sep="\t", header=None)
        counts = df[2].value_counts()

        for alg, n in counts.items():
            rows.append({
                "m": m,
                "ALG": alg,
                "markers": n
            })

    if not rows:
        raise RuntimeError("No input files found")

    df = pd.DataFrame(rows)
    df = df.sort_values(["m", "markers"], ascending=[True, False])

    # --------------------------------------------------
    # Plot
    # --------------------------------------------------

    color_cycle = cycle(OKABE_ITO)

    fig, ax = plt.subplots(figsize=(16, 7))
    bottom = defaultdict(int)

    for alg in df["ALG"].unique():
        sub = df[df["ALG"] == alg]
        ms = sub["m"]
        heights = sub["markers"]
        bottoms = [bottom[m] for m in ms]

        ax.bar(ms, heights, bottom=bottoms, color=next(color_cycle), label=str(alg))

        for m, h in zip(ms, heights):
            bottom[m] += h

    ax.set_xlabel("m")
    ax.set_ylabel(f"Number of markers at {node}")
    ax.set_title(f"{prefix}: marker composition of {node} across m")

    ax.set_xlim(df["m"].min() - 0.5, df["m"].max() + 0.5)

    plt.tight_layout()
    plt.savefig(f"{out}.pdf")
    plt.savefig(f"{out}.png", dpi=300)

    print(f"Saved {out}.pdf")
    print(f"Saved {out}.png")


if __name__ == "__main__":
    main()
