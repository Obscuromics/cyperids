#!/usr/bin/env python3

import os
import csv
import glob

CHUNK_SIZE = 5000

# You are running this from inside CAP_run
BASE_DIR = "."

CHUNKS_DIR = os.path.join(BASE_DIR, "chunks_nuclear") #output from satellite_windows.py
HELIXER_DIR = os.path.join(BASE_DIR, "Helixer_annotations")
RM_DIR = os.path.join(BASE_DIR, "RepeatMasker_outputs")
OUT_DIR = os.path.join(BASE_DIR, "chunks_nuclear_augmented_CDS")

os.makedirs(OUT_DIR, exist_ok=True)


def chunk_coords(chunk_id):
    start = (chunk_id - 1) * CHUNK_SIZE + 1
    end = chunk_id * CHUNK_SIZE
    return start, end


def overlaps(a_start, a_end, b_start, b_end):
    return not (a_end < b_start or b_end < a_start)


def read_helixer_cds(gff_path):
    cds = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, feature = parts[0], parts[2]
            if feature != "CDS":
                continue
            start, end = int(parts[3]), int(parts[4])
            cds.setdefault(chrom, []).append((start, end))
    return cds


def is_TE(class_family):
    non_te = {"Simple_repeat", "Low_complexity", "Satellite", "Unknown"}
    return class_family not in non_te


def read_repeatmasker_out(rm_path):
    repeats = {}
    with open(rm_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("SW") or line.startswith("score"):
                continue
            parts = line.split()
            chrom = parts[4]
            start = int(parts[5])
            end = int(parts[6])
            class_family = parts[10]
            if not is_TE(class_family):
                continue
            repeats.setdefault(chrom, []).append((start, end, class_family))
    return repeats


chunk_files = glob.glob(os.path.join(CHUNKS_DIR, "*_chunks.csv"))
if not chunk_files:
    raise RuntimeError("No chunk files found in chunks_nuclear/")

for chunk_file in chunk_files:
    species = os.path.basename(chunk_file).replace("_chunks.csv", "")
    print(f"Processing {species}")

    helixer_gff = os.path.join(HELIXER_DIR, f"{species}.helixer.gff3")
    rm_out = os.path.join(RM_DIR, f"{species}.fna.out")

    if not os.path.exists(helixer_gff):
        print(f"  Helixer file missing for {species}, skipping")
        continue
    if not os.path.exists(rm_out):
        print(f"  RepeatMasker file missing for {species}, skipping")
        continue

    cds = read_helixer_cds(helixer_gff)
    repeats = read_repeatmasker_out(rm_out)

    out_name = f"{species}_chunks_augmented_CDS.csv"
    out_path = os.path.join(OUT_DIR, out_name)

    with open(chunk_file) as infile, open(out_path, "w", newline="") as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ["CDS_overlap", "TE_overlap", "TE_class"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            chrom = row["chromosome"]
            chunk_id = int(row["chunk_id"])
            c_start, c_end = chunk_coords(chunk_id)

            cds_overlap = 0
            for cds_start, cds_end in cds.get(chrom, []):
                if overlaps(c_start, c_end, cds_start, cds_end):
                    cds_overlap = 1
                    break

            te_classes = set()
            for r_start, r_end, cls in repeats.get(chrom, []):
                if overlaps(c_start, c_end, r_start, r_end):
                    te_classes.add(cls)

            row["CDS_overlap"] = cds_overlap
            row["TE_overlap"] = 1 if te_classes else 0
            row["TE_class"] = ";".join(sorted(te_classes)) if te_classes else ""

            writer.writerow(row)

print("Done.")
