#!/bin/bash
# Write cell volumes and cell centres into each family's control case at time 0.
# flow_structure.py reads V/Cx/Cy/Cz from there; the mesh is shared within a
# family (verified by md5 of constant/polyMesh/points), so once per family is enough.
set +eu
source /usr/lib/openfoam/openfoam2506/etc/bashrc >/dev/null 2>&1
set -u

for f in "$HOME/results_matched_0.2T_f4000hz" \
         "$HOME/openfoam_projects/results" \
         "$HOME/openfoam_projects/results_spot" \
         "$HOME/openfoam_projects/results_weld_s5e6"; do
    [ -d "$f/no-mhd" ] || { echo "MISSING $f/no-mhd"; continue; }
    cd "$f/no-mhd" || continue
    postProcess -func writeCellVolumes -time 0 >/dev/null 2>&1
    postProcess -func writeCellCentres -time 0 >/dev/null 2>&1
    n=0
    for g in V Cx Cy Cz; do [ -f "0/$g" ] && n=$((n + 1)); done
    echo "$n/4 geometry fields  $f"
done
