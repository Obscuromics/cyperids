#!/usr/bin/env bash
set -euo pipefail

#############################################
# USAGE
#############################################

if [[ $# -lt 4 ]]; then
    echo "Usage:"
    echo "  $0 PREFIX NODE ALG_PLOT_SCRIPT STACK_PLOT_SCRIPT INPUT_DIR M_START [M_END]"
fi

PREFIX="$1"
NODE="$2"
ALG_PLOT_SCRIPT="$(realpath "$3")"
STACK_PLOT_SCRIPT="$(realpath "$4")"
INPUT_DIR="$5"
M_START="$6"
M_END="${7:-$6}"   # if M_END not given, do just one m

#############################################
# MAIN LOOP
#############################################

for m in $(seq "${M_START}" "${M_END}"); do
    echo "Processing m=${m}"

    workdir="m${m}_output"
    tab="tabulate_m${m}"

    mkdir -p "${workdir}"

    # copy inputs if missing
    if ! ls "${workdir}/${PREFIX}${m}"* >/dev/null 2>&1; then
        cp "${INPUT_DIR}/${PREFIX}${m}"* "${workdir}/"
    fi

    cd "${workdir}"

    PICKLE="${PREFIX}${m}.with_ancestors.pickle"
    if [[ ! -f "${PICKLE}" ]]; then
        echo "WARNING: missing ${PICKLE}, skipping m=${m}"
        cd ..
        continue
    fi

    syngraph tabulate \
        -g "${PICKLE}" \
        -o "${tab}"

    # -----------------------------------------
    # Extract BUSCO-style table for given node
    # -----------------------------------------

    awk -F'\t' -v node="${NODE}_" '
    NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i ~ "^" node) cols[i]=1
        }
    }
    NR>1 {
        for (i in cols) {
            if ($i != "NA") {
                print $1, "Complete", $i, 0, 1
            }
        }
    }
    ' OFS='\t' "${tab}.table.tsv" \
    > "m${m}.buscona_${NODE}.${PREFIX}.tsv"

    # -----------------------------------------
    # Per-m ALG number plot
    # -----------------------------------------

    python3 "${ALG_PLOT_SCRIPT}" \
        "m${m}.buscona_${NODE}.${PREFIX}.tsv" \
        -o "m${m}_${NODE}_${PREFIX}_markercounts"

    cd ..

done


    #############################################
    # FINAL STACKED PLOT ACROSS ALL m
    #############################################

python3 "${STACK_PLOT_SCRIPT}" \
    --prefix "${PREFIX}" \
    --node "${NODE}" \
    --out "${PREFIX}_${NODE}_evaluation"
