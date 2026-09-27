# Experimental pore births — ESRF, AlSi10Mg 0.5 mm

In-situ X-ray radiography, recorded 2025-11-14. **6,898 tracked pore births**,
6,578 inside the laser window. Produced by `github.com/Teslos/yolo-pores`
(`export_pore_births.py -d results/tiff`); the working copy lives in
`mhd-openfoam/experiment/` and these are curated copies.

| file | what |
|---|---|
| `pore_births.csv` | one row per pore |
| `pore_weights.csv` | normalised radial / depth / diameter histograms in SI metres, per condition and pooled as `ALL` |
| `pore_births.md` | column dictionary and the full caveat list — **read it before fitting anything** |

Nine recordings: 100 W / 100 ms **spot** welds, 50,400 fps, 2.2 µm/px. One
`nomf` pair plus one each of `bx bxy bxz by byz bz`, all 400 Hz alternating.
Chosen because it is the only block with a verified 0.00 false-positive floor
over 13,655 pre-laser frames.

## The axis convention is a cyclic permutation — check `sim_label`

The `condition` letters are the **experimental** convention. Simulation labels
differ, and the two planes that contain the beam axis swap:

| physical | experiment | simulation |
|---|---|---|
| transverse | `bx` | B∥z |
| scan | `by` | B∥x |
| vertical | `bz` | B∥y |
| horizontal plane | `bxy` | **XZ** |
| transverse-vertical plane | `bxz` | **YZ** |
| scan-vertical plane | `byz` | **XY** |

Both CSVs carry `field_direction` and `sim_label` columns. Use those, never the
letters.

## What has been used so far

`cases/poreCloud-rmf-xy-0.2T-expbyz` takes its `radialSeeding` from condition
`byz` (= sim XY, 689 births). `scripts/test_poreSeedingFromData.py` regenerates
those weights from this directory, so the dictionary and the data cannot drift
apart unnoticed.

**Only the radial distribution was transferred.** The diameters were not, and
should not be without a deliberate decision: mean measured pore volume is
3.257e-13 m³ against 5.59e-15 m³ in the matched SS316L VoF runs — **58× larger**
— which with the calibrated injection rate implies a 252% void fraction against
the 0.83–2.53% the resolved runs measure. Aluminium at 100 W spot welds simply
makes bigger pores than steel at 150 W scanning, and since the exclusion force
scales as d³, importing the sizes would confound alloy with field.

## Caveats that travel with any use of this data

* **Regime.** Spot welds, aluminium. The simulations are 150 W at 36 mm/s on
  SS316L. There is no scan, so a radius here is about a *static* axis.
* **`birth` is first detection, not nucleation.** Recall ≈ 0.55 and a pore is
  invisible inside the keyhole, so births are systematically late and the radial
  histogram is biased **outward**.
* **No z.** Radiography projects through the plate, so `r_um` is the in-plane
  radius and a lower bound — tight, since 0.5 mm is the dimension along the beam.
* **Diameters are truncated with pile-up.** Detector floor 20.6 µm; the 25–30 µm
  bin holds 23% of all pores. Do not fit below ~35–40 µm.
* **Time origin is good to ~20 ms.** Prefer `t_from_kh_end_ms` over `birth_t_ms`.
* **`--min-life 5`**, not the 20 of the published fate tables; 43.7% of rows
  would survive that cut. Filter on `life_frames` to match.

## Regenerating a seeding block

    ~/.venvs/mhd-paper-analysis/bin/python3 scripts/poreSeedingFromData.py \
        <a CSV of r_um,d_um for one condition> --col r=r_um --scale-um \
        --edges 0 25e-6 50e-6 75e-6 100e-6 125e-6 150e-6 175e-6 200e-6 225e-6 250e-6 \
        --i-have-birth-data

Pick the outer edge from the melt pool, not the data: the XY pool reaches 230 µm
at the t = 1 ms restart, so 250 µm wastes no weight. Bins beyond the liquid get
`binMult = 0` in `MeltPoolInjection.C` and their weight is silently discarded,
lowering the effective injection rate.

## Not included

The 1 mm campaign — 6 mm/s scanning welds on a re-used plate carrying 2–7
pre-existing pores per frame, so most of its "births" are old pores drifting
into view.
