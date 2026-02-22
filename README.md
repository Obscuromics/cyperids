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
