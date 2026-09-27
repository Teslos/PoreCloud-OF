#!/bin/bash
# Flow-structure scoring across every case family that has a matched no-field control.
#
#   results_matched_0.2T_f4000hz  the richest set: 3 DC axes + 3 RMF planes + control,
#                                 all at 0.2 T / 4000 Hz. PRIMARY.
#   results                       0.1 and 0.2 T of each plane -> carries the chaos bound
#   results_spot                  spot-weld geometry, same 0.1/0.2 T set
#   results_weld_s5e6             the high-sigma (5 MS/m) variant
#
# Result families are not in git (they are ~440 GB). Point RESULTS_ROOT and
# MATCHED_ROOT at wherever they live on this machine; the defaults are the
# layout on the machine this study was run on.
set -u

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_ROOT="${RESULTS_ROOT:-$HOME/openfoam_projects}"
MATCHED_ROOT="${MATCHED_ROOT:-$HOME/results_matched_0.2T_f4000hz}"
out="${OUT:-$HOME/results/flow_structure_report.txt}"

mkdir -p "$(dirname "$out")"
: > "$out"

run() {
    echo "##############################################################" >> "$out"
    echo "########## $1" >> "$out"
    echo "##############################################################" >> "$out"
    shift
    if [ ! -d "$2" ]; then
        echo "SKIP: $2 not present on this machine" >> "$out"; echo >> "$out"; return
    fi
    python3 "$here/flowStructureSurvey.py" "$@" >> "$out" 2>&1
    echo >> "$out"
}

run "PRIMARY: 0.2 T matched, 4000 Hz -- 3 DC axes + 3 RMF planes vs no field" \
    --family "$MATCHED_ROOT" --max-snaps 30

run "0.1/0.2 T sweep (weld) -- includes the CHAOS BOUND" \
    --family "$RESULTS_ROOT/results" --max-snaps 30 --chaos

run "0.1/0.2 T sweep (spot geometry) -- includes the CHAOS BOUND" \
    --family "$RESULTS_ROOT/results_spot" --max-snaps 30 --chaos

run "high-sigma 5 MS/m variant -- includes the CHAOS BOUND" \
    --family "$RESULTS_ROOT/results_weld_s5e6" --max-snaps 30 --chaos

echo "DONE -> $out"
wc -l "$out"
