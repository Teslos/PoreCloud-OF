# Experimental pore data — drop files here

Nothing is here yet. The seeding currently shipped in `cases/*/constant/poreCloudProperties`
is derived from the **resolved-VoF runs**, not from experiment — see
`../poretracker-summaries/`. This directory is where experimental data replaces it.

## What to put here

A CSV of pore **births**, one row per pore. Any column names; you map them on
the command line. The useful columns:

| quantity | needed for |
|---|---|
| birth position — either a radius from the laser axis, or x and z | `radialSeeding.weights` |
| birth time | only if the laser scans, so the axis position is known per pore |
| diameter | `sizeDistribution` |
| depth below the free surface | not consumed yet; record it if you have it |

Then:

    python3 ../../scripts/poreSeedingFromData.py yourfile.csv \
        --col x=X_um --col z=Z_um --col t=t_ms --col d=diam_um --scale-um \
        --melt-volume 1.70e-11 --window 1e-3 \
        --i-have-birth-data

It prints the `rate`, `radialSeeding` and `sizeDistribution` blocks ready to
paste into `constant/poreCloudProperties`. `scripts/test_poreSeedingFromData.py`
anchors it by regenerating the weights currently in use from the VoF summaries,
so the script and the dictionaries cannot drift apart unnoticed.

## The one distinction that decides whether the data is usable

**Micro-CT of a finished track is the wrong measurement.** It records where the
*surviving* pores came to rest — births already filtered by flotation, transport
and capture at the solidification front. That filtering is exactly what a
simulation seeded from this data would be trying to predict, so using it is
circular and biased toward survivors.

**In-situ X-ray radiography is the right one.** It times and locates nucleation
events, which is what the injector needs.

The script cannot tell them apart, so it refuses to guess: pass
`--i-have-birth-data` or `--final-positions`. The second still works, and stamps
the warning into the emitted comment so the provenance travels with the numbers.

## Match the regime, not just the alloy

These runs are **150 W at 36 mm/s** — slow, keyhole-welding conditions. Typical
LPBF radiography is taken at 500–1000 mm/s, where pore formation is a different
regime. Matching power and scan speed matters more than matching the alloy.

## Why this would be worth doing

Three known weaknesses in the current VoF-derived seeding that experimental data
would fix:

1. **The size distribution has a detection floor.** PoreTracker's
   `MIN_PORE_CELLS 5` piles pores up at ~9.3 µm — which is why p5, p10 and p25
   of that data are identical — and anything below ~12 µm is invisible to it.
2. **The radial distribution is strongly field-dependent** and therefore not
   transferable between configurations: χ² = 111.9 on 30 df, p < 0.0001 across
   the seven matched cases. The weights in use pool all seven, which is
   defensible for `rmf-yz` (χ² = 8.1, p = 0.15 against the others) and **not**
   for the rest — `rmf-xz` puts 37.9% of births in the 80–100 µm bin against the
   pooled 19.3%.
3. **It inherits whatever the VoF solver gets wrong** about where pores form.

Point 2 is why per-case experimental data, or per-case PoreTracker passes, are
needed before seeding any configuration other than `rmf-yz`.
