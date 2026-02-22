# cyperids

#### fission_fusion_counter_plus.sh
This script takes _summary.tsv files from https://github.com/charlottewright/lep_busco_painter.
Firstly, BUSCO tables for internal nodes need to be fabricated, having the same format as the BUSCO output full_table.tsv files for a typical genome. The position coordinates are unimportant, and chromosome name should be substituted for the ALG. Internal nodes can now be used as either "reference" or "query" genomes for lep_busco_painter. The _summary.tsv files should be generated for every pair of "reference" and "query" genomes, where the "reference" is the parent node and the "query" is the child node.

#### clade_rearrangement_comparisons.R
This script takes rearrangement numbers and a newick tree (see /data) to compare rearrangement rates between Cyperaceae and Juncaceae, and between Carex and non-Carex Cyperaceae.

#### heterozygosity.R
This script compares heterozygosity to rearrangement rates along terminal branches (input files in /data)

#### modelling_fusion_probability.R
This script runs the simulations/MLE for the four models described in Supplementary Text 2, and includes the plotting scripts for the relevant visualisations. Input files are found in /data.

#### round_1-4_consensus_building.sh
Concatenates RepeatModeler round 1-4 consensi.

#### satellite_jaccards.R
Computes 8-mer Jaccard similarity between consensus sequences for each candidate oligocentromere family per species (see /data for consensus sequences).

#### satellite_windows.py
Uses CAP results, a csv file of chromosome lengths, and a text file with the CAP repeat class names of the candidate centromeric satellites (see /data) in order to score, for each 5 kb window along the chromosome, the overlap of a candidate centromeric satellite array.

#### TE_CDS_windows.py
Creates augmented versions of the csv files produced with satellite_windows.py to include CDS and TE overlap from Helixer and RepeatModeler annotations.

#### oligocentromere_enrichment_and_decay.R
Takes scores of CDS overlap, TE overlap, and candidate centromeric satellite overlap of 5 kb windows from each species (see /data/window-data) and calculates and visualises TE enrichment and gene decay with proximity to the centromeres.


