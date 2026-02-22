import csv
from pathlib import Path
from collections import defaultdict

# ----------------------------
# Paths and parameters
# ----------------------------
BASE = Path(".")
CAP_RESULTS = BASE / "CAP_results"
CHROM_LEN = BASE / "chromosome_lengths"
OUTDIR = BASE / "chunks"
OUTDIR.mkdir(exist_ok=True)

CHUNK_SIZE = 5000

# ----------------------------
# Read centromeric repeats table
# ----------------------------
centro_repeats = {}

with open(BASE / "satellite_repeats.txt", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        species = row["Species"].replace(" ", "_")
        reps = row["Centromeric repeat"].replace('"', "")
        centro_repeats[species] = {r.strip() for r in reps.split(",")}

# ----------------------------
# Process each species
# ----------------------------
for species_dir in CAP_RESULTS.iterdir():
    if not species_dir.is_dir():
        continue

    species = species_dir.name
    if species not in centro_repeats:
        continue

    print(f"Processing {species}")

    repeats_of_interest = centro_repeats[species]

    # ----------------------------
    # Read chromosome lengths
    # ----------------------------
    chrom_lengths = {}
    chrom_len_file = CHROM_LEN / f"{species}.chromosome_lengths.csv"

    if not chrom_len_file.exists():
        print(f"  WARNING: missing chromosome lengths for {species}")
        continue

    with open(chrom_len_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom_lengths[row["chromosome"]] = int(float(row["length_bp"]))

    # ----------------------------
    # Read repeat instances
    # ----------------------------
    repeats = defaultdict(list)
    repeats_file = species_dir / f"{species}.fna_repeats_with_seq_filtered_reclassed.csv"

    if not repeats_file.exists():
        print(f"  WARNING: missing repeats file for {species}")
        continue

    with open(repeats_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rep_class = row["new_class"]
            if rep_class not in repeats_of_interest:
                continue

            try:
                start = int(float(row["start"]))
                end = int(float(row["end"]))
            except ValueError:
                continue

            chrom = row["seqID"].replace('"', "")
            repeats[chrom].append((start, end, rep_class))

    # ----------------------------
    # Write chunked output
    # ----------------------------
    out_file = OUTDIR / f"{species}_chunks.csv"

    with open(out_file, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["chromosome", "chunk_id", "repeat_classes"])

        for chrom, length in chrom_lengths.items():
            n_chunks = (length + CHUNK_SIZE - 1) // CHUNK_SIZE

            for i in range(n_chunks):
                chunk_start = i * CHUNK_SIZE + 1
                chunk_end = min((i + 1) * CHUNK_SIZE, length)

                hits = set()
                for r_start, r_end, rep in repeats.get(chrom, []):
                    if r_end >= chunk_start and r_start <= chunk_end:
                        hits.add(rep)

                writer.writerow([
                    chrom,
                    i + 1,
                    ";".join(sorted(hits)) if hits else ""
                ])

    print(f"  → wrote {out_file}")
