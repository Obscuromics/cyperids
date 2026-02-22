#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage:"
    echo "  $0 <input_tsv_dir> <tree.newick> <output.csv>"
    exit 1
fi

INDIR="$1"
TREE="$2"
OUTCSV="$3"

# Step 1: Compute fissions/fusions per TSV

tmp_counts=$(mktemp)

for f in "$INDIR"/*.tsv; do
    fname=$(basename "$f")

    parent="${fname%%_*}"
    child="${fname#*_}"
    child="${child%_summary.tsv}"

    read -r fissions fusions < <(
        awk '
        {
            row_gt = 0
            for (i=1; i<=NF; i++) {
                if ($i > 10) {
                    row_gt++
                    col_gt[i]++
                }
            }
            if (row_gt > 0) {
                fusions += (row_gt - 1)
            }
        }
        END {
            for (i in col_gt) {
                if (col_gt[i] > 0) {
                    fissions += (col_gt[i] - 1)
                }
            }
            print fissions+0, fusions+0
        }' "$f"
    )

    rearrangements=$((fissions + fusions))

    echo "$parent,$child,$fissions,$fusions,$rearrangements" >> "$tmp_counts"
done

# Step 2: Attach branch lengths + compute rates

python3 <<EOF
from ete3 import Tree
import csv

tree = Tree("$TREE", format=1)

# Build lookup: (parent, child) -> branch_length
branch = {}
for node in tree.traverse():
    if node.up and node.name and node.up.name:
        branch[(node.up.name, node.name)] = node.dist

with open("$OUTCSV", "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow([
        "parent",
        "child",
        "fissions",
        "fusions",
        "rearrangements",
        "branch_length",
        "fission_rate",
        "fusion_rate",
        "rearrangement_rate"
    ])

    with open("$tmp_counts") as inp:
        for line in inp:
            parent, child, fissions, fusions, rearr = line.strip().split(",")

            fissions = float(fissions)
            fusions = float(fusions)
            rearr = float(rearr)

            bl = branch.get((parent, child), None)

            if bl is None or bl == 0:
                fission_rate = ""
                fusion_rate = ""
                rearr_rate = ""
                bl_out = ""
            else:
                fission_rate = fissions / bl
                fusion_rate = fusions / bl
                rearr_rate = rearr / bl
                bl_out = bl

            writer.writerow([
                parent,
                child,
                int(fissions),
                int(fusions),
                int(rearr),
                bl_out,
                fission_rate,
                fusion_rate,
                rearr_rate
            ])
EOF

rm "$tmp_counts"

echo "Summary written to $OUTCSV"
