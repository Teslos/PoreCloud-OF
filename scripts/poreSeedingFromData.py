#!/usr/bin/env python3
"""Turn a table of pore births into the poreCloudProperties seeding block.

Three things in `injectionModels/melt` describe where pores come from, and all
three are currently derived from the resolved-VoF runs rather than from
experiment:

    rate              births per unit LIQUID volume per second
    radialSeeding     birth radius about the (moving) laser axis
    sizeDistribution  pore diameters

This script computes all three from one table of pore births and prints the
dictionary text, so replacing the VoF proxy with experimental data is a paste
rather than a hand calculation.

    python3 poreSeedingFromData.py births.csv --melt-volume 1.70e-11

Column names are given with --col, so the script does not care whose export it
is reading:

    python3 poreSeedingFromData.py exp.csv \\
        --col r=radius_um --col d=diameter_um --col t=time_ms --scale-um

WHAT THE INPUT MUST BE
----------------------
Pore *births*: where and when each pore first appeared in the liquid.

Micro-CT of a finished track is NOT this. It records the final resting place of
the pores that survived, which is births already filtered by flotation,
transport and capture at the solidification front - the very thing a simulation
using this input is trying to predict. Feeding it in is circular, and biased
toward whatever survives. In-situ radiography, which times and locates
nucleation events, is the measurement that fits.

The script cannot detect the difference, so it asks: pass --i-have-birth-data
to confirm, or --final-positions to proceed anyway with the caveat recorded in
the emitted comment.
"""

import argparse
import csv
import math
import sys

import numpy as np

# Default radial bins, in metres. Matches the bins the VoF-derived weights use,
# so a new dataset is directly comparable to data/poretracker-summaries/.
DEFAULT_EDGES = [0, 2e-5, 4e-5, 6e-5, 8e-5, 1e-4, 1.3e-4]

# The diameter table is emitted on a UNIFORM grid because OpenFOAM's `general`
# distribution with `cumulative false` integrates it TRAPEZOIDALLY as a density
# (general.C:57-63). Supplying probability MASS on unevenly spaced points
# over-weights every wide bin by its width - doing exactly that once put the
# mean bubble volume 2.9x too high.
DIAM_STEP = 2e-6


def load(path, cols, scale):
    """Read the table, returning columns as float arrays under canonical names."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{path}: no data rows")

    have = set(rows[0])
    out = {}
    for key, name in cols.items():
        if name is None:
            continue
        if name not in have:
            raise SystemExit(
                f"{path}: no column {name!r} (for --col {key}=). "
                f"Available: {', '.join(sorted(have))}"
            )
        vals = []
        for r in rows:
            v = r[name]
            vals.append(float(v) if v not in ("", None) else math.nan)
        out[key] = np.asarray(vals) * scale
    return out, len(rows)


def radius(d, axis_xz, axis_v):
    """Birth radius about the laser axis, which moves with the scan."""
    if "r" in d:
        return d["r"]
    for k in ("x", "z"):
        if k not in d:
            raise SystemExit(
                "need either --col r=<radius column>, or both --col x= and "
                "--col z= so the radius can be computed"
            )
    x0, z0 = axis_xz
    t = d.get("t")
    if t is None:
        if axis_v:
            raise SystemExit(
                "--axis-velocity given but no --col t=<birth time>; the axis "
                "moves, so the radius is undefined without it"
            )
        t = np.zeros_like(d["x"])
    return np.hypot(d["x"] - (x0 + axis_v * t), d["z"] - z0)


def weights(r, edges):
    r = r[np.isfinite(r)]
    inside = r < edges[-1]
    cnt, _ = np.histogram(r[inside], bins=edges)
    if cnt.sum() == 0:
        raise SystemExit(f"no births inside {edges[-1]*1e6:.0f} um")
    return cnt, cnt / cnt.sum(), int((~inside).sum())


def diameters(dia):
    """Histogram diameters onto the uniform grid `general` needs."""
    dia = dia[np.isfinite(dia) & (dia > 0)]
    if dia.size == 0:
        return None
    lo = math.floor(dia.min() / DIAM_STEP) * DIAM_STEP
    hi = math.ceil(dia.max() / DIAM_STEP) * DIAM_STEP + DIAM_STEP
    edges = np.arange(lo, hi + DIAM_STEP, DIAM_STEP)
    cnt, _ = np.histogram(dia, bins=edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    dens = cnt / (cnt.sum() * DIAM_STEP)          # a density, per the note above

    # Every bin is emitted, empty ones included. Dropping zeros looks tidier and
    # is wrong: the abscissae would no longer be uniformly spaced, so the
    # trapezoidal integration would run a straight line ACROSS each gap instead
    # of through zero, quietly restoring the very width-weighting the uniform
    # grid exists to prevent.
    return centres, dens * DIAM_STEP, dia


def emit(src, n_in, cnt, w, dropped, edges, dia, rate, provenance, caveat):
    p = print
    p("// " + "-" * 68)
    p(f"// Seeding derived from {src}")
    p(f"// {provenance}")
    if caveat:
        for line in caveat:
            p(f"// {line}")
    p(f"// {n_in} rows in, {cnt.sum()} births inside {edges[-1]*1e6:.0f} um"
      + (f", {dropped} dropped beyond it" if dropped else ""))
    p("// " + "-" * 68)
    if rate is not None:
        p(f"rate            {rate:.6g};")
    p("radialSeeding")
    p("{")
    p("    axisPoint     (250e-6 0 250e-6);")
    p("    axisVelocity  (0.036 0 0);")
    p("    edges   ( " + " ".join(f"{e:g}" for e in edges) + " );")
    p("    weights ( " + " ".join(f"{x:.5f}" for x in w) + " );")
    p("}")
    p(f"// per-bin counts: {cnt.tolist()}")

    if dia is None:
        p("// no diameter column given - sizeDistribution left unchanged")
        return
    centres, mass, raw = dia
    p("")
    p("// Densities on a UNIFORM 2 um grid, NOT per-bin probabilities:")
    p("// `general` with cumulative false integrates trapezoidally as a")
    p("// density (general.C:57-63), so mass on uneven spacing over-weights")
    p("// wide bins. Uniform spacing removes the ambiguity.")
    p(f"// n = {raw.size}, mean d = {raw.mean()*1e6:.2f} um, "
      f"mean volume = {(math.pi/6*raw**3).mean():.3e} m3")
    p("sizeDistribution")
    p("{")
    p("    type            general;")
    p("    generalDistribution")
    p("    {")
    p("        cumulative      false;")
    p("        distribution")
    p("        (")
    for c, m in zip(centres, mass):
        p(f"            ({c:.3e}   {m:.5f})")
    p("        );")
    p("    }")
    p("}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", help="CSV of pore births")
    ap.add_argument("--col", action="append", default=[], metavar="KEY=NAME",
                    help="map a canonical key (r, x, z, y, d, t) to a column "
                         "name. Repeatable.")
    ap.add_argument("--scale-um", action="store_true",
                    help="input lengths and times are in um / us, not m / s")
    ap.add_argument("--edges", type=float, nargs="+", default=DEFAULT_EDGES,
                    metavar="M", help="radial bin edges in metres")
    ap.add_argument("--axis-point", type=float, nargs=2, default=(250e-6, 250e-6),
                    metavar=("X", "Z"), help="laser axis at t=0, metres")
    ap.add_argument("--axis-velocity", type=float, default=0.036,
                    metavar="MPS", help="scan speed along x")
    ap.add_argument("--melt-volume", type=float, default=None, metavar="M3",
                    help="time-averaged liquid volume; with --window gives rate")
    ap.add_argument("--window", type=float, default=None, metavar="S",
                    help="duration the births were counted over")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--i-have-birth-data", action="store_true",
                   help="confirm the table records nucleation, not final state")
    g.add_argument("--final-positions", action="store_true",
                   help="proceed with final/CT positions, recording the caveat")
    args = ap.parse_args(argv)

    if not (args.i_have_birth_data or args.final_positions):
        ap.error(
            "say which this is: --i-have-birth-data (in-situ radiography, "
            "nucleation events) or --final-positions (micro-CT of a finished "
            "track). They are not interchangeable - see the module docstring.")

    cols = {"r": None, "x": None, "z": None, "y": None, "d": None, "t": None}
    for spec in args.col:
        if "=" not in spec:
            ap.error(f"--col wants KEY=NAME, got {spec!r}")
        k, v = spec.split("=", 1)
        if k not in cols:
            ap.error(f"unknown --col key {k!r}; pick from {', '.join(cols)}")
        cols[k] = v
    if not any(cols[k] for k in ("r", "x")):
        ap.error("need --col r= or --col x= and --col z=")

    scale = 1e-6 if args.scale_um else 1.0
    d, n_in = load(args.table, cols, scale)
    r = radius(d, args.axis_point, args.axis_velocity)
    cnt, w, dropped = weights(r, args.edges)
    dia = diameters(d["d"]) if "d" in d else None

    rate = None
    if args.melt_volume and args.window:
        rate = cnt.sum() / (args.melt_volume * args.window)

    caveat = []
    if args.final_positions:
        caveat = [
            "WARNING: built from FINAL (e.g. micro-CT) positions, not births.",
            "These are births already filtered by transport and capture, so",
            "seeding with them is circular and biased toward survivors.",
        ]
    prov = ("in-situ birth events" if args.i_have_birth_data
            else "final positions - see warning")
    emit(args.table, n_in, cnt, w, dropped, args.edges, dia, rate, prov, caveat)
    return 0


if __name__ == "__main__":
    sys.exit(main())
