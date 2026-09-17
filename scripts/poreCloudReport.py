#!/usr/bin/env python3
"""Report on a poreCloud run, reusing mhd-openfoam's plotting layer.

The heavy part - seven pages of matplotlib, the summary table, the publication
palette - already exists in mhd-openfoam/analyse_pore_forces.py, written for
pores found by connected components in the resolved-VoF fields.  This script
does not reimplement any of it.  It only converts poreCloud_pores.csv into the
track dicts that make_pdf() already consumes, and calls it.

Two things make the conversion much simpler than the VoF path:

  - Parcels carry a persistent PoreID, so grouping replaces the greedy
    nearest-neighbour matching that _build_pore_tracks() needs.
  - Forces are reported per bubble, so no field sampling is needed.

Unit contract, checked in test_poreCloudReport.py:

  CSV                          track dict                 note
  Cx_m, Cy_m, Cz_m  [m]        x, depth, width  [um]      *1e6; 'depth' is the
                                                          y coordinate, matching
                                                          analyse_pore_forces
  Time  [s]                    t_us  [us]                 *1e6
  Fx_N..  [N]                  EF  [N/m3]                 / Volume_m3
  Fbuoyx.. [N]                 LF  [N/m3]                 / Volume_m3
  Ux.. - Uslip..  [m/s]        U_pore  [m/s]              carrier velocity:
                                                          Uslip = U_parcel - Uc
  Volume_m3                    size  [cells]              / cell_volume

The force slots are the report's, not ours: 'EF' holds the Leenov-Kolin
exclusion force, which is what the VoF exclusionForce field holds too, so that
slot means the same thing in both. 'LF' holds buoyancy rather than the melt's
Lorentz force, which the CSV does not carry - hence the label override.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# Default mesh of the 0.2 T azimuthal case: 500x700x500 um over 60x96x48 cells.
DEFAULT_CELL_VOLUME = (500e-6 / 60) * (700e-6 / 96) * (500e-6 / 48)

# Serial coverage floor for that case (docs/specs.md 5.5).  Resolution
# dependent: take it from a serial run before trusting it elsewhere.
DEFAULT_COVERAGE_FLOOR = 0.89


def _import_report_layer(explicit=None):
    """Import mhd-openfoam's analyse_pore_forces, or explain how to point at it."""
    candidates = [
        explicit,
        os.environ.get('MHD_OPENFOAM_DIR'),
        Path.home() / 'mhd-openfoam',
    ]
    for c in candidates:
        if c and (Path(c) / 'analyse_pore_forces.py').exists():
            sys.path.insert(0, str(Path(c)))
            import analyse_pore_forces as apf
            return apf
    raise SystemExit(
        'Cannot find analyse_pore_forces.py (the plotting layer this reuses).\n'
        'Pass --report-lib /path/to/mhd-openfoam or set MHD_OPENFOAM_DIR.'
    )


def load_rows(csv_path):
    """Read the CSV, newest schema only (a Coverage column must be present)."""
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f'{csv_path} has no data rows')
    if 'Coverage' not in rows[0]:
        raise SystemExit(
            f'{csv_path} predates the Coverage column - rerun on a library '
            'build from 9adc766 or later.'
        )
    return rows


def build_tracks(rows, cell_volume=DEFAULT_CELL_VOLUME, min_len=4,
                 active_only=True):
    """Group CSV rows into analyse_pore_forces-style track dicts.

    active_only drops samples of parcels already entrapped by the
    solidification front: they are held at a fixed position by design, so
    including them would report a velocity of zero as if it were physics.
    """
    by_id = defaultdict(list)
    for r in rows:
        if active_only and int(float(r['Active'])) != 1:
            continue
        by_id[r['PoreID']].append(r)

    tracks = []
    for pid, recs in by_id.items():
        recs.sort(key=lambda r: float(r['Time']))
        if len(recs) < min_len:
            continue

        vol = np.array([float(r['Volume_m3']) for r in recs])
        vol_safe = np.where(vol > 0, vol, np.nan)

        def vec(*cols):
            return np.array([[float(r[c]) for c in cols] for r in recs])

        td = {
            'pore_id': pid,
            't_us':   np.array([float(r['Time']) for r in recs]) * 1e6,
            'x':      np.array([float(r['Cx_m']) for r in recs]) * 1e6,
            'depth':  np.array([float(r['Cy_m']) for r in recs]) * 1e6,
            'width':  np.array([float(r['Cz_m']) for r in recs]) * 1e6,
            'size':   vol / cell_volume,
            'T_pore': np.array([float(r['T_K']) for r in recs]),
            'coverage': np.array([float(r['Coverage']) for r in recs]),
            'd_over_dx': np.array([float(r['dOverDx']) for r in recs]),
            # Force densities, so the units match the report's N/m3 axes.
            'EF': vec('Fx_N', 'Fy_N', 'Fz_N') / vol_safe[:, None],
            'LF': vec('Fbuoyx', 'Fbuoyy', 'Fbuoyz') / vol_safe[:, None],
            # Carrier velocity: Uslip = U_parcel - Uc  =>  Uc = U - Uslip
            'U_pore': vec('Ux', 'Uy', 'Uz') - vec('Uslipx', 'Uslipy', 'Uslipz'),
        }

        _add_kinematics(td)
        tracks.append(td)

    tracks.sort(key=lambda t: t['t_us'][0])
    return tracks


def _add_kinematics(td):
    """Velocities, net displacement, lifetime and impulses.

    Central differences interior, one-sided at the ends - the same scheme
    analyse_pore_forces.analyse_case uses, so the two reports stay comparable.
    Positions are um and times us, so velocity is um/us == m/s.
    """
    t = td['t_us']
    n = len(t)
    v = {}
    for key, axis in (('vx', 'x'), ('vdepth', 'depth'), ('vwidth', 'width')):
        a = td[axis]
        d = np.zeros(n)
        if n >= 3:
            dt = t[2:] - t[:-2]
            d[1:-1] = np.where(dt > 0, (a[2:] - a[:-2]) / np.where(dt > 0, dt, 1), 0)
        if n >= 2:
            d[0] = (a[1] - a[0]) / (t[1] - t[0]) if t[1] > t[0] else 0.0
            d[-1] = (a[-1] - a[-2]) / (t[-1] - t[-2]) if t[-1] > t[-2] else 0.0
        v[key] = d
    td.update(v)
    td['vmag'] = np.sqrt(v['vx'] ** 2 + v['vdepth'] ** 2 + v['vwidth'] ** 2)

    td['displacement'] = np.array([
        td['x'][-1] - td['x'][0],
        td['depth'][-1] - td['depth'][0],
        td['width'][-1] - td['width'][0],
    ])
    td['lifetime_us'] = t[-1] - t[0]

    dt_s = np.diff(t) * 1e-6
    td['LF_impulse'] = np.nansum(td['LF'][:-1] * dt_s[:, None], axis=0)
    td['EF_impulse'] = np.nansum(td['EF'][:-1] * dt_s[:, None], axis=0)
    td['total_impulse'] = td['LF_impulse'] + td['EF_impulse']


def coverage_stats(rows, floor=DEFAULT_COVERAGE_FLOOR):
    """Truncation statistics - the quantity specific to a decomposed run."""
    cov = np.array([float(r['Coverage']) for r in rows])
    return {
        'n': len(cov),
        'mean': float(cov.mean()),
        'median': float(np.median(cov)),
        'min': float(cov.min()),
        'frac_below_floor': float((cov < floor).mean()),
        'frac_below_half': float((cov < 0.5).mean()),
    }


def capture_stats(rows):
    """Final-snapshot split of the cloud into moving bubbles and frozen pores."""
    t_end = max(float(r['Time']) for r in rows)
    last = [r for r in rows if abs(float(r['Time']) - t_end) < 1e-15]
    active = sum(1 for r in last if int(float(r['Active'])) == 1)
    return {
        't_end': t_end,
        'n_parcels': len(last),
        'n_active': active,
        'n_pores': len(last) - active,
    }


def force_discontinuity(tracks, jump=0.10, steady=0.01, dx=None):
    """Is a coverage change accompanied by a step in force?

    Returns per-displacement-bin medians of |dF|/F, split by whether the
    stencil changed between the two samples.  Binning by distance travelled is
    the control: a bubble whose stencil changed also moved further, and would
    show a larger force change for that reason alone.
    """
    dx = dx or DEFAULT_CELL_VOLUME ** (1 / 3)
    edges = [(0, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 1.0), (1.0, np.inf)]
    bins = {e: ([], []) for e in edges}

    for td in tracks:
        f = np.linalg.norm(td['EF'], axis=1)
        pos = np.stack([td['x'], td['depth'], td['width']], axis=1) * 1e-6
        step = np.linalg.norm(np.diff(pos, axis=0), axis=1) / dx
        dcov = np.abs(np.diff(td['coverage']))
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = np.abs(np.diff(f)) / f[:-1]

        for i in range(len(rel)):
            if not np.isfinite(rel[i]) or f[i] <= 0:
                continue
            for e in edges:
                if e[0] <= step[i] < e[1]:
                    if dcov[i] < steady:
                        bins[e][0].append(rel[i])
                    elif dcov[i] > jump:
                        bins[e][1].append(rel[i])
                    break

    out = []
    for e in edges:
        s, j = bins[e]
        if len(s) < 50 or len(j) < 50:
            continue
        ms, mj = float(np.median(s)), float(np.median(j))
        out.append({
            'bin': e, 'n_steady': len(s), 'n_jump': len(j),
            'median_steady': ms, 'median_jump': mj,
            'ratio': mj / ms if ms > 0 else float('nan'),
        })
    return out


def print_report(rows, tracks, disc):
    cs, ps = coverage_stats(rows), capture_stats(rows)
    print(f"\nrows {len(rows):,}   tracks {len(tracks):,}")
    print(f"final cloud at t={ps['t_end']:.6f}s: {ps['n_parcels']} parcels = "
          f"{ps['n_active']} active + {ps['n_pores']} pores "
          f"({ps['n_pores'] / max(ps['n_parcels'], 1):.0%} captured)")
    print(f"\ncoverage  mean {cs['mean']:.3f}  median {cs['median']:.3f}  "
          f"min {cs['min']:.3f}")
    print(f"          below floor {cs['frac_below_floor']:.1%}   "
          f"below 0.50 {cs['frac_below_half']:.1%}")

    if disc:
        print(f"\nforce step at stencil changes, controlled for displacement:")
        print(f"{'displacement':>16} {'steady':>9} {'jump':>9} {'ratio':>7}")
        for d in disc:
            lo, hi = d['bin']
            name = f'{lo:g}-{hi:g} dx' if np.isfinite(hi) else f'>{lo:g} dx'
            print(f"{name:>16} {d['median_steady'] * 100:8.2f}% "
                  f"{d['median_jump'] * 100:8.2f}% {d['ratio']:6.1f}x")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv', type=Path, help='poreCloud_pores.csv')
    ap.add_argument('-o', '--output', type=Path, help='PDF path')
    ap.add_argument('--report-lib', help='directory holding analyse_pore_forces.py')
    ap.add_argument('--label', default='poreCloud', help='case label in the legend')
    ap.add_argument('--min-len', type=int, default=4, help='min samples per track')
    ap.add_argument('--max-tracks', type=int, default=40,
                    help='longest-lived N tracks to plot (0 = all); the report '
                         'draws one line per track and is unreadable past a few dozen')
    ap.add_argument('--cell-volume', type=float, default=DEFAULT_CELL_VOLUME)
    ap.add_argument('--coverage-floor', type=float, default=DEFAULT_COVERAGE_FLOOR)
    args = ap.parse_args()

    apf = _import_report_layer(args.report_lib)
    # The 'LF' slot holds buoyancy here, not the melt's Lorentz force.
    apf.LF_LABEL, apf.LF_TITLE = 'Buoyancy', 'Buoyancy force at bubble'
    apf.EF_LABEL, apf.EF_TITLE = 'Exclusion', 'Leenov-Kolin force at bubble'

    rows = load_rows(args.csv)
    tracks = build_tracks(rows, args.cell_volume, args.min_len)
    if not tracks:
        raise SystemExit(f'no track reached {args.min_len} samples')

    disc = force_discontinuity(tracks)
    print_report(rows, tracks, disc)

    plotted = sorted(tracks, key=lambda t: -t['lifetime_us'])
    if args.max_tracks:
        plotted = plotted[:args.max_tracks]

    apf.print_summary_table([plotted], [args.label])
    out = args.output or args.csv.with_suffix('.pdf')
    apf.make_pdf([plotted], [args.label], out)
    print(f'\nreport -> {out}  ({len(plotted)} of {len(tracks)} tracks plotted)')


if __name__ == '__main__':
    main()
