#!/bin/bash
# Provision a fresh workstation for the MHD / pore-expulsion modelling stack.
#
# Code only: every repo is cloned, every binary is rebuilt. None of the ~440 GB
# of simulation output comes along. The case dictionaries and calibration inputs
# live in this repo (cases/, data/), so nothing has to be copied from the old
# machine -- this script needs only git access.
#
# Run it from a clone of PoreCloud-OF:
#
#   git clone git@github.com:Teslos/PoreCloud-OF.git
#   ./PoreCloud-OF/scripts/bootstrapWorkstation.sh
#
#   ... verify     just re-run the checks
#
# Deliberately NOT set -e: the apt and build phases each report their own
# status, and one missing optional package should not abort the rest.
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELFREPO="$(cd "$HERE/.." && pwd)"

# This script ships inside PoreCloud-OF, so it is normally run from a clone that
# already exists. Adopt that clone as the canonical one instead of cloning a
# second copy into $PROJ and then symlinking cases/ and data/ at the wrong one.
if [ -d "$SELFREPO/.git" ] && [ -f "$SELFREPO/src/poreCloud/poreCloudMagneticField.C" ]; then
    PCREPO="$SELFREPO"
else
    PCREPO=""          # set to $PROJ/PoreCloud-OF once provision_repos clones it
fi
FOAM_BASHRC=/usr/lib/openfoam/openfoam2506/etc/bashrc
PROJ="$HOME/openfoam_projects"
VENV="$HOME/.venvs/mhd-paper-analysis"
fail=0

# OpenFOAM's bashrc reads unset variables, so `set -u` turns sourcing it into a
# silent abort. Always go through this wrapper.
foam_env() {
    set +u
    # shellcheck disable=SC1090
    source "$FOAM_BASHRC"
    local rc=$?
    set -u
    return $rc
}

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }
ok()  { printf '  ok    %s\n' "$*"; }
bad() { printf '  FAIL  %s\n' "$*"; fail=$((fail+1)); }

# ---------------------------------------------------------------- 1. OpenFOAM
provision_openfoam() {
    say "OpenFOAM v2506 (apt, dl.openfoam.com)"
    if [ -f "$FOAM_BASHRC" ]; then ok "already installed"; return; fi
    curl -fsSL https://dl.openfoam.com/add-debian-repo.sh | sudo bash || {
        bad "could not add the OpenFOAM apt repo -- check network/proxy"; return; }
    sudo apt-get update
    sudo apt-get install -y openfoam2506-default openfoam2506-dev \
                            openfoam2506-source openfoam2506-tools
    [ -f "$FOAM_BASHRC" ] && ok "installed" || bad "install did not produce $FOAM_BASHRC"
}

# ------------------------------------------------------------- 2. system deps
provision_apt() {
    say "build toolchain and python/plotting packages"
    # ffmpeg was MISSING on the old machine, which made the animation steps in
    # visualize_*.py fail silently. Install it here.
    sudo apt-get install -y build-essential openmpi-bin libopenmpi-dev \
        python3-venv python3-numpy python3-scipy python3-matplotlib \
        imagemagick ffmpeg git curl rsync
}

# -------------------------------------------------------------- 3. shell env
provision_shellrc() {
    say "shell environment"
    # The old machine's rc also sourced /tmp/rc_fhn_uv_bin/env -- a stale shim
    # that does not survive a reboot. It is deliberately not reproduced.
    if grep -qF "$FOAM_BASHRC" "$HOME/.bashrc" 2>/dev/null; then
        ok ".bashrc already sources OpenFOAM"
    else
        printf '\n# OpenFOAM v2506\nsource %s\n' "$FOAM_BASHRC" >> "$HOME/.bashrc"
        ok "added the OpenFOAM source line to .bashrc"
    fi
}

# --------------------------------------------------------------- 4. venv
provision_venv() {
    say "analysis venv ($VENV)"
    # Rebuilt, never copied: a venv bakes absolute paths into its scripts.
    # --system-site-packages is load-bearing -- numpy/scipy/matplotlib come from
    # the system, and only these two are venv-local. System python3 has NO
    # pandas, so anything importing pandas must run under this interpreter.
    [ -d "$VENV" ] || python3 -m venv --system-site-packages "$VENV"
    "$VENV/bin/pip" install -q --upgrade pip
    "$VENV/bin/pip" install -q pandas==3.0.3 pymupdf==1.28.0 && ok "pandas + pymupdf"
}

# --------------------------------------------------------------- 5. clones
clone() {  # clone <url> <dest> <branch>
    if [ -d "$3/.git" ]; then
        ok "$(basename "$3") present; pulling"
        git -C "$3" pull --ff-only 2>&1 | tail -1
    else
        # --depth 1: the solver's history alone is 295 MB (it contains a 95 MB
        # committed log file) against ~36 MB of tracked source. This box is for
        # computing; run `git fetch --unshallow` if it becomes a dev machine.
        git clone --depth 1 --branch "$2" "$1" "$3" || bad "clone failed: $1"
    fi
}

provision_repos() {
    say "repositories"
    mkdir -p "$PROJ"
    clone git@github.com:Coolnesss/mhd-laserbeamfoam-solver.git azimuthal-field "$HOME/mhd-laserbeamfoam-solver"
    if [ -n "$PCREPO" ]; then
        ok "PoreCloud-OF: using the clone this script was run from ($PCREPO)"
        git -C "$PCREPO" pull --ff-only 2>&1 | tail -1
    else
        PCREPO="$PROJ/PoreCloud-OF"
        clone git@github.com:Teslos/PoreCloud-OF.git master "$PCREPO"
    fi
    clone git@github.com:Teslos/PoreTracker-OF.git             main            "$PROJ/PoreTracker-OF"

    # mhd-openfoam tracks a 46 MB constant/polyMesh that blockMesh regenerates,
    # so it is excluded at checkout rather than downloaded and deleted.
    local m="$HOME/mhd-openfoam"
    if [ ! -d "$m/.git" ]; then
        git clone --depth 1 --branch experimental --filter=blob:none --sparse \
            git@ssh.gitlab.empa.ch:intelligent-manufacturing-group/mhd-openfoam.git "$m" \
          && git -C "$m" sparse-checkout set --no-cone '/*' '!/constant/polyMesh' \
          || bad "mhd-openfoam clone failed (Empa GitLab reachable? SSH key present?)"
    else
        ok "mhd-openfoam present"
    fi
}

# --------------------------------------------------------------- 6. build
build() {
    say "build"
    foam_env || { bad "cannot source $FOAM_BASHRC"; return; }
    ( cd "$HOME/mhd-laserbeamfoam-solver" && ./Allwmake -j ) 2>&1 | tail -3
    ( cd "$PCREPO"                        && ./Allwmake -j ) 2>&1 | tail -3
    # PoreTracker-OF has no Allwmake.
    ( cd "$PROJ/PoreTracker-OF"           && wmake libso src/PoreTracker ) 2>&1 | tail -3
}

# --------------------------------------------------------------- 7. inputs
install_inputs() {
    say "calibration inputs and case dictionaries"
    # These live in the repo, so there is nothing to transfer. Symlinked rather
    # than copied: a copy would drift from the clone the moment either changed,
    # and these are inputs under version control, not scratch data.
    mkdir -p "$HOME/results"
    ln -sfn "$PCREPO/data"  "$HOME/results/poreCloud-data"
    ln -sfn "$PCREPO/cases" "$HOME/results/poreCloud-cases"
    ok "~/results/poreCloud-{data,cases} -> the clone"
    printf '  note  cases carry system/ + constant/ only. To run one:\n'
    printf '          cd <case> && blockMesh && cp -r ../initial-common 0 && setFields\n'
    printf '  note  mhd-openfoam working-tree changes are NOT in git. If the old\n'
    printf '        machine had any, move them across yourself.\n'
}

# --------------------------------------------------------------- 8. verify
verify() {
    say "verification"
    foam_env || true

    [ "${WM_PROJECT_VERSION:-}" = v2506 ] \
        && ok "WM_PROJECT_VERSION=v2506" || bad "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"

    for x in laserbeamFoam poreCloudReplay setSolidFraction; do
        [ -x "$FOAM_USER_APPBIN/$x" ] && ok "app $x" || bad "app $x missing"
    done
    for l in libPoreCloud.so libPoreTracker.so liblaserHeatSource.so; do
        [ -f "$FOAM_USER_LIBBIN/$l" ] && ok "lib $l" || bad "lib $l missing"
    done

    # These two catch the one failure mode a green build would hide: a clone
    # taken before the old machine pushed. Both features were committed last,
    # so their absence means the source is stale even though everything builds.
    grep -q alternating_MHD_value \
        "$HOME/mhd-laserbeamfoam-solver/applications/solvers/laserbeamFoam/UEqn.H" \
        && ok "AC field mode present" || bad "AC field mode MISSING -- stale clone"
    grep -q radialSeeding \
        "$PCREPO/src/poreCloud/injection/MeltPoolInjection.C" \
        && ok "measured-distribution seeding present" || bad "radialSeeding MISSING -- stale clone"

    # The cloud builds its own B. Before this branch existed it silently used a
    # rotating field for an AC case, giving a spurious vertical exclusion force
    # where the applied field produces exactly none (specs 5.13). A clone from
    # before the fix builds and runs cleanly, so only the source shows it.
    grep -q "alternating_" "$PCREPO/src/poreCloud/poreCloudMagneticField.C" \
        && ok "cloud honours MHD.alternating" \
        || bad "cloud AC branch MISSING -- stale clone, AC runs will be wrong"

    ( cd "$PCREPO" && python3 scripts/test_poreCloudReport.py >/dev/null 2>&1 ) \
        && ok "poreCloudReport tests" || bad "poreCloudReport tests"
    ( cd "$PCREPO" && python3 scripts/test_poreCloudFieldSurvey.py >/dev/null 2>&1 ) \
        && ok "poreCloudFieldSurvey tests" || bad "poreCloudFieldSurvey tests"

    "$VENV/bin/python3" -c 'import pandas, fitz' 2>/dev/null \
        && ok "venv: pandas + pymupdf import" || bad "venv imports"

    command -v ffmpeg >/dev/null && ok "ffmpeg" || bad "ffmpeg (animations will fail silently)"

    printf '\n'
    [ "$fail" -eq 0 ] && printf '\033[1mall checks passed\033[0m\n' \
                      || printf '\033[1m%d check(s) FAILED\033[0m\n' "$fail"
    return "$fail"
}

if [ "${1:-all}" = verify ]; then
    [ -n "$PCREPO" ] || PCREPO="$PROJ/PoreCloud-OF"
    verify; exit $?
fi
provision_openfoam
provision_apt
provision_shellrc
provision_venv
provision_repos
build
install_inputs
verify
