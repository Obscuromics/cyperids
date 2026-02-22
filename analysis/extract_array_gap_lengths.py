#!/usr/bin/env python3

import csv
import glob
import os

CHUNK_SIZE_KB = 5

def process_file(filepath):
    species = os.path.basename(filepath).replace("_chunks.csv", "")
    results = []

    with open(filepath, newline="") as f:
        reader = csv.DictReader(f)
        rows = sorted(reader, key=lambda r: (r["chromosome"], int(r["chunk_id"])))

    current_chr = None
    current_type = None
    run_length = 0

    def flush_run():
        if current_chr is None:
            return
        results.append({
            "species": species,
            "chromosome": current_chr,
            "type": current_type,
            "chunks": run_length,
            "length_kb": run_length * CHUNK_SIZE_KB
        })

    for row in rows:
        chrom = row["chromosome"]
        is_array = row["repeat_classes"].strip() != ""
        run_type = "array" if is_array else "gap"

        if chrom != current_chr:
            flush_run()
            current_chr = chrom
            current_type = run_type
            run_length = 1
        elif run_type != current_type:
            flush_run()
            current_type = run_type
            run_length = 1
        else:
            run_length += 1

    flush_run()
    return results


def main():
    files = glob.glob("*_chunks.csv")
    if not files:
        raise RuntimeError("No *_chunks.csv files found")

    all_results = []

    for file in files:
        all_results.extend(process_file(file))

    with open("array_gap_lengths.csv", "w", newline="") as out:
        fieldnames = ["species", "chromosome", "type", "chunks", "length_kb"]
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)


if __name__ == "__main__":
    main()
