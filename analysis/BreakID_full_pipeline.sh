#!/bin/bash
# BreakID full pipeline with optional precomputed minimap2 alignments
# Usage:
# ./BreakID_full_pipeline.sh genome1.fa genome2.fa busco1.tsv busco2.tsv output_prefix [--paf1 paf1] [--paf2 paf2]

set -e
set -o pipefail

GENOME1=$1
GENOME2=$2
BUSCO1=$3
BUSCO2=$4
OUTPREFIX=$5
shift 5

# Optional arguments
PAF1=""
PAF2=""
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --paf1) PAF1="$2"; shift 2 ;;
    --paf2) PAF2="$2"; shift 2 ;;
    *) shift ;;
  esac
done

WIN_SIZE=50000
STEP_SIZE=10000
MEDIAN_SIZE=101

OUTDIR_CAREX1=${OUTPREFIX}_$GENOME1
OUTDIR_CAREX2=${OUTPREFIX}_$GENOME2
mkdir -p $OUTDIR_CAREX1
mkdir -p $OUTDIR_CAREX2

echo "=============================="
echo "Step 0: Minimap2 alignments"
echo "=============================="

if [[ -z "$PAF1" || -z "$PAF2" ]]; then
    echo "No precomputed PAFs supplied; running minimap2..."
    PAF1=${OUTDIR_CAREX1}/genome2_vs_genome1.paf
    PAF2=${OUTDIR_CAREX2}/genome1_vs_genome2.paf

    echo "Aligning Genome1 -> Genome2..."
    minimap2 -cx asm10 $GENOME1 $GENOME2 > $PAF1
    echo "Aligning Genome2 -> Genome1..."
    minimap2 -cx asm10 $GENOME2 $GENOME1 > $PAF2
else
    echo "Using provided PAF files:"
    echo "PAF1: $PAF1"
    echo "PAF2: $PAF2"
fi

echo "=============================="
echo "Step 1: BreakID Part 1 (coarse breakpoints)"
echo "=============================="

run_breakid_part1() {
    local PAF=$1
    local OUTDIR=$2
    mkdir -p $OUTDIR
    echo "Processing $PAF -> $OUTDIR"

    # Chromosome sizes
    cut -f1-2 $PAF | sort | uniq > $OUTDIR/genome_query.genome
    cut -f6-7 $PAF | sort | uniq > $OUTDIR/genome_target.genome

    # Extract alignment positions
    cut -f1,3,4,6,8,9 $PAF | sort -k4,4 -k5,5n > $OUTDIR/aln_pairs.txt

    # Assign numeric IDs to query chromosomes
    cut -f1 $OUTDIR/genome_query.genome | nl -n ln > $OUTDIR/query_chr_num.txt
    awk 'NR==FNR {map[$2]=$1; next} { $1=map[$1]; print $4"\t"$5"\t"$6"\t"$1}' \
        $OUTDIR/query_chr_num.txt $OUTDIR/aln_pairs.txt > $OUTDIR/target_genome_composition.bed

    # Sliding windows
    bedtools makewindows -g $OUTDIR/genome_target.genome -w $WIN_SIZE -s $STEP_SIZE > $OUTDIR/windows.bed
    bedtools map -a $OUTDIR/windows.bed -b $OUTDIR/target_genome_composition.bed -c 4 -o mean -null 0 \
        > $OUTDIR/target_genome_composition_smoothed.bed

    # Coarse breakpoints: robust R code
    Rscript - <<EOF
library(dplyr)
library(stats)

MEDIAN_SIZE <- $MEDIAN_SIZE
df <- read.table("$OUTDIR/target_genome_composition_smoothed.bed", header=FALSE)
colnames(df) <- c("chr","start","end","chrnum")

breaks <- data.frame()
for(chr in unique(df\$chr)){
    df_chr <- df[df\$chr==chr & df\$chrnum != 0,]
    if(nrow(df_chr)==0) next
    df_chr\$medianSig <- runmed(df_chr\$chrnum, MEDIAN_SIZE, endrule="constant")
    df_chr\$change_to_prev <- c(0, diff(df_chr\$medianSig))
    potential_cp <- df_chr[df_chr\$change_to_prev != 0,]

    if(nrow(potential_cp) == 0){
        next
    } else if(nrow(potential_cp) == 1){
        breaks <- rbind(breaks, data.frame(chr=chr, group=nrow(breaks)+1,
                                           break_start=potential_cp\$start[1],
                                           break_end=potential_cp\$end[1]))
        next
    }

    interval_s <- potential_cp\$start[1]
    interval_e <- potential_cp\$end[1]

    for(i in 2:nrow(potential_cp)){
        if(is.na(potential_cp\$start[i]) | is.na(potential_cp\$end[i-1])) next
        if(potential_cp\$start[i] - potential_cp\$end[i-1] < 1e6){
            interval_e <- potential_cp\$end[i]
        } else {
            breaks <- rbind(breaks, data.frame(chr=chr, group=nrow(breaks)+1,
                                               break_start=interval_s,
                                               break_end=interval_e))
            interval_s <- potential_cp\$start[i]
            interval_e <- potential_cp\$end[i]
        }
    }
    breaks <- rbind(breaks, data.frame(chr=chr, group=nrow(breaks)+1,
                                       break_start=interval_s,
                                       break_end=interval_e))
}

write.table(breaks, file="$OUTDIR/break_intervals.txt", sep="\t", row.names=FALSE, quote=FALSE)
EOF
}

run_breakid_part1 $PAF1 $OUTDIR_CAREX1
run_breakid_part1 $PAF2 $OUTDIR_CAREX2

echo "=============================="
echo "Step 2: Classify breakpoints + robust PDF"
echo "=============================="

classify_breaks() {
    local BREAK_FILE=$1
    local BUSCO_FILE=$2
    local OUTDIR=$3
    local TARGET_SPECIES=$4
    local QUERY_SPECIES=$5

    Rscript - <<EOF
library(dplyr)

# --- FIXED: robust BUSCO reading with proper column specification ---
busco <- read.table("$BUSCO_FILE", header=TRUE, stringsAsFactors=FALSE, 
                    sep="\t", colClasses=c("character", "character", "numeric", "character", "character"))

# Ensure proper column names
colnames(busco) <- c("buscoID", "query_chr", "position", "assigned_chr", "status")

# Debug: print structure
cat("BUSCO file structure:\n")
str(busco)
cat("\nFirst few rows:\n")
print(head(busco))

# Verify position is numeric
if(!is.numeric(busco\$position)) {
    stop("Position column is not numeric after reading!")
}

# Read breaks file
breaks <- read.table("$BREAK_FILE", header=TRUE, stringsAsFactors=FALSE, sep="\t")
breaks\$break_start <- as.numeric(breaks\$break_start)
breaks\$break_end <- as.numeric(breaks\$break_end)

# compute BUSCO color changes
busco_changes <- busco %>%
  arrange(query_chr, position) %>%
  group_by(query_chr) %>%
  mutate(prev_assigned = lag(assigned_chr, default=first(assigned_chr)),
         color_change = assigned_chr != prev_assigned) %>%
  filter(color_change) %>%
  select(query_chr, position) %>%
  ungroup() %>%
  as.data.frame()

if(nrow(busco_changes)==0) busco_changes <- data.frame(query_chr=character(0), position=numeric(0))

cat("\nNumber of BUSCO color changes detected:", nrow(busco_changes), "\n")

classify_break <- function(chr, start, end, busco_changes_chr, chr_in_busco){
  # If chromosome not in BUSCO data at all
  if(!chr_in_busco) return("undefined")
  
  # If chromosome is in BUSCO but has no color changes
  if(nrow(busco_changes_chr)==0) return(paste("fission in", "$TARGET_SPECIES"))
  
  pos <- as.numeric(busco_changes_chr\$position)
  dist_to_change <- abs(pos - ((start+end)/2))
  min_dist <- if(length(dist_to_change) > 0) min(dist_to_change, na.rm=TRUE) else Inf
  
  # Close to BUSCO color change = fusion in focal/target species
  # Far from BUSCO color change = fission in focal/target species  
  if(min_dist <= 5e5){
    return(paste("fusion in", "$QUERY_SPECIES"))
  } else {
    return(paste("fission in", "$TARGET_SPECIES"))
  }
}

breaks\$class <- sapply(1:nrow(breaks), function(i){
    chr <- breaks\$chr[i]
    start <- breaks\$break_start[i]
    end <- breaks\$break_end[i]
    chr_in_busco <- chr %in% unique(busco\$query_chr)
    busco_chr <- busco_changes %>% filter(query_chr==chr)
    classify_break(chr, start, end, busco_chr, chr_in_busco)
})

OUTFILE="${OUTDIR}/break_intervals_classified.txt"
write.table(breaks, file=OUTFILE, sep="\t", row.names=FALSE, quote=FALSE)

cat("\nClassification complete. Summary:\n")
print(table(breaks\$class))

# PDF: one chromosome per page
pdf(sub("break_intervals_classified.txt","breaks_classified_signal.pdf", OUTFILE),
    width=10, height=4)
chrlist <- unique(breaks\$chr)
for(chr in chrlist){
    df_chr <- breaks[breaks\$chr==chr,]
    if(nrow(df_chr)==0) next
    x <- (df_chr\$break_start + df_chr\$break_end)/2
    y <- 1:nrow(df_chr)
    plot(x, y, type="n", xlab="Chromosome position", ylab="", main=chr)
    cols <- ifelse(grepl("fusion", df_chr\$class), "red",
                   ifelse(grepl("fission", df_chr\$class), "blue", "gray"))
    points(x, y, col=cols, pch=19)
    legend("topright", legend=c("fusion","fission","undefined"), col=c("red","blue","gray"), pch=19)
}
dev.off()

cat("\nPDF created successfully\n")
EOF
}

# classify breaks for both genomes
classify_breaks ${OUTDIR_CAREX1}/break_intervals.txt $BUSCO1 $OUTDIR_CAREX1 $GENOME1 $GENOME2
classify_breaks ${OUTDIR_CAREX2}/break_intervals.txt $BUSCO2 $OUTDIR_CAREX2 $GENOME2 $GENOME1

echo "=============================="
echo "Done!"
echo "Outputs:"
echo "${OUTDIR_CAREX1}/break_intervals_classified.txt"
echo "${OUTDIR_CAREX1}/breaks_classified_signal.pdf"
echo "${OUTDIR_CAREX2}/break_intervals_classified.txt"
echo "${OUTDIR_CAREX2}/breaks_classified_signal.pdf"
echo "=============================="
