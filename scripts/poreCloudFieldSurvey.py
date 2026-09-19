#!/usr/bin/env python3
"""Where does the exclusion force point? An Eulerian survey of stored fields.

The Lagrangian route (run the cloud, track bubbles, integrate forces) costs
hours of CFD per case and, on a finished run, cannot be replayed: the velocity
field decorrelates ~60% between stored snapshots, so interpolating between them
reconstructs nothing. See docs/specs.md.

But the question that matters - does the Leenov-Kolin exclusion force push
bubbles toward the free surface? - is a property of the INSTANTANEOUS field,
not of trajectories. It needs no time integration and therefore no
interpolation, so it can be answered exactly from any finished run.

The physics that makes this cheap:

    F_exclusion = -1.5 * V_bubble * (J x B) = -1.5 * V_bubble * LorentzForce

so the exclusion force is ANTIPARALLEL to the melt's own Lorentz body force,
which laserbeamFoam already writes. Verified on stored fields: the solver's
own exclusionForce sits at cos(angle) = -0.9993 to LorentzForce. No B field,
no cloud, no interpolation - just a sign flip on a field already on disk.

Validation: on azimuthal-0.2T this gives 30.3% of the liquid pushed upward,
against 29.6% of bubbles from the full Lagrangian run (specs 5.8). Two
independent methods, one percentage point apart.

What it does NOT give: trajectories, residence times, or whether a bubble
actually escapes. Those need the Lagrangian cloud on finely-written fields.

Usage:
    python3 poreCloudFieldSurvey.py <case> [<case> ...] [-o out.csv] [--stride N]
"""

import argparse
import csv
import os
import re
import sys

import numpy as np

# Free surface: metal occupies y below this in these cases. Expulsion is +y.
FREE_SURFACE_M = 547e-6
ALPHA_MIN = 0.5
EPS_MIN = 0.5


def read_field(path, vec):
    """Read an OpenFOAM internalField, ascii or binary, uniform or nonuniform.

    Format is taken from the FoamFile header, NOT guessed. An earlier version
    tried the binary interpretation whenever the byte count happened to fit and
    accepted it if the values looked finite and not absurd. That guard is far
    weaker than it appears: ascii digit bytes reinterpreted as little-endian
    float64 usually land on small, finite, entirely plausible numbers. It
    silently returned 1e-259..1e-33 for a temperature field whose real range is
    306..3892 K, which made a `T > 870` mask empty and an analysis built on it
    wrong. Scalars are the dangerous case - a vector needs three times the
    bytes, so ascii is too short to be mistaken for one.
    """
    raw = open(path, 'rb').read()
    kind = b'vector' if vec else b'scalar'

    fmt = re.search(rb'^\s*format\s+(ascii|binary)\s*;', raw[:2048], re.M)
    binary = bool(fmt) and fmt.group(1) == b'binary'

    m = re.search(rb'internalField\s+nonuniform\s+List<' + kind + rb'>\s*\n(\d+)\s*\n\(', raw)
    if m:
        n = int(m.group(1))
        start = m.end()
        count = n * (3 if vec else 1)
        if binary:
            a = np.frombuffer(raw[start:start + 8 * count], dtype='<f8')
            if len(a) != count:
                raise ValueError(f'{path}: truncated binary block')
            return a.reshape(n, 3) if vec else a
        if raw[start:start + 1] == b'\n':
            start += 1
        end = raw.index(b'\n)', start)
        toks = raw[start:end].replace(b'(', b' ').replace(b')', b' ').split()
        a = np.array([float(x) for x in toks])
        return a.reshape(n, 3) if vec else a

    m2 = re.search(rb'internalField\s+uniform\s+\(?([-0-9.eE+ ]+?)\)?\s*;', raw)
    if m2:
        v = [float(x) for x in m2.group(1).split()]
        if vec:
            return np.tile(v, (1, 1))          # caller broadcasts
        return np.array(v)
    raise ValueError(f'cannot parse internalField in {path}')


def stored_times(case):
    """Time directory names that actually hold the fields we need, in order."""
    out = []
    for d in os.listdir(case):
        if not re.fullmatch(r'[0-9.eE+-]+', d):
            continue
        if all(os.path.exists(os.path.join(case, d, f))
               for f in ('LorentzForce', 'alpha.metal', 'epsilon1')):
            out.append(d)
    return sorted(out, key=float)


def survey_snapshot(case, t, ncells=None):
    """One snapshot: how much of the liquid has the exclusion force pointing up?

    Returns None when the pool is too small to mean anything - a nearly
    solidified or not-yet-formed pool gives a percentage computed over a
    handful of cells, which is noise wearing a number's clothes.
    """
    d = os.path.join(case, t)
    alpha = read_field(os.path.join(d, 'alpha.metal'), False)
    eps = read_field(os.path.join(d, 'epsilon1'), False)
    lf = read_field(os.path.join(d, 'LorentzForce'), True)

    n = max(len(alpha), len(eps), len(lf))
    if ncells:
        n = ncells
    if len(alpha) == 1:
        alpha = np.full(n, alpha[0])
    if len(eps) == 1:
        eps = np.full(n, eps[0])
    if lf.shape[0] == 1:
        lf = np.tile(lf[0], (n, 1))

    liquid = (alpha > ALPHA_MIN) & (eps > EPS_MIN)
    if liquid.sum() < 100:
        return None

    # Exclusion force is antiparallel to the Lorentz body force.
    ef_y = -lf[liquid][:, 1]
    mag = np.linalg.norm(lf[liquid], axis=1)
    up = ef_y > 0

    # Two ways the question can have no answer, both of which must report NaN
    # rather than 0%, because 0% ranks a case as the WORST configuration when
    # it is really no configuration at all.
    #
    #   1. No field: F = 0 everywhere, so there is no direction to ask about.
    #
    #   2. B parallel to y: the Lorentz force is J x B, which is perpendicular
    #      to B by construction, so F_y is IDENTICALLY zero - not small, exactly
    #      zero in every cell - while F_x and F_z are large. Such a field cannot
    #      push a bubble up or down even in principle. Measured on
    #      matched/dc-by-0.2T: F_y == 0 in all 57306 cells with force, against
    #      max|F_z| = 1.28e8. Reporting that as "0% pushed up" invites the
    #      reading "pushes everything down", which is the opposite of true:
    #      it exerts no vertical force at all.
    degenerate_y = mag.max() > 0 and np.all(ef_y == 0.0)
    if mag.max() <= 0 or degenerate_y:
        return {
            'time': float(t), 'n_liquid': int(liquid.sum()),
            'frac_up': float('nan'), 'frac_up_force_weighted': float('nan'),
            'median_force_density': float(np.median(mag)) if mag.max() > 0 else 0.0,
            'median_ef_y': 0.0,
            'undefined_reason': 'B_parallel_to_y' if degenerate_y else 'no_field',
        }

    return {
        'time': float(t),
        'n_liquid': int(liquid.sum()),
        'frac_up': float(up.mean()),
        'frac_up_force_weighted':
            float((mag * up).sum() / mag.sum()) if mag.sum() > 0 else float('nan'),
        'median_force_density': float(np.median(mag)),
        'median_ef_y': float(np.median(ef_y)),
        'undefined_reason': '',
    }


def survey_case(case, stride=10):
    times = stored_times(case)
    if not times:
        return None, 'no time directories with LorentzForce/alpha.metal/epsilon1'
    picks = times[::max(1, len(times) // stride)] or times
    if times[-1] not in picks:
        picks.append(times[-1])

    rows = []
    for t in picks:
        try:
            r = survey_snapshot(case, t)
        except Exception as e:                       # noqa: BLE001
            return None, f'{type(e).__name__} at t={t}: {e}'
        if r:
            rows.append(r)
    if not rows:
        return None, 'no snapshot had a resolvable melt pool'
    return rows, None


def summarise(case, rows):
    fu = np.array([r['frac_up'] for r in rows])
    fw = np.array([r['frac_up_force_weighted'] for r in rows])
    w = np.array([r['n_liquid'] for r in rows], float)     # weight by pool size
    nofield = bool(np.isnan(fu).all())
    return {
        'case': case,
        'n_snapshots': len(rows),
        'no_field': nofield,
        'undefined_reason': next((r.get('undefined_reason','') for r in rows
                                  if r.get('undefined_reason')), ''),
        'frac_up_mean': float('nan') if nofield else float(np.nanmean(fu)),
        'frac_up_poolweighted': float('nan') if nofield
            else float(np.nansum(fu * w) / np.nansum(w * ~np.isnan(fu))),
        'frac_up_force_weighted': float(np.nanmean(fw)),
        'frac_up_final': float(fu[-1]),
        'median_force_density': float(np.median([r['median_force_density'] for r in rows])),
        'max_liquid_cells': int(w.max()),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cases', nargs='+')
    ap.add_argument('-o', '--output', default=None, help='write per-case summary CSV')
    ap.add_argument('--per-snapshot', default=None, help='write every snapshot row here')
    ap.add_argument('--stride', type=int, default=10, help='snapshots per case (default 10)')
    args = ap.parse_args()

    summaries, snaps, failures = [], [], []
    for case in args.cases:
        label = case.rstrip('/').replace('/home/ivt/', '')
        rows, err = survey_case(case.rstrip('/'), args.stride)
        if err:
            failures.append((label, err))
            print(f'SKIP {label}: {err}', file=sys.stderr)
            continue
        s = summarise(label, rows)
        summaries.append(s)
        for r in rows:
            snaps.append(dict(case=label, **r))

    summaries.sort(key=lambda s: -s['frac_up_poolweighted'])
    print(f'\n{"#":>3} {"case":<46} {"snaps":>6} {"up(pool-wtd)":>13} '
          f'{"up(force-wtd)":>14} {"|F| N/m3":>11}')
    for i, s in enumerate(summaries, 1):
        print(f'{i:3d} {s["case"]:<46} {s["n_snapshots"]:6d} '
              f'{s["frac_up_poolweighted"]:12.1%} {s["frac_up_force_weighted"]:13.1%} '
              f'{s["median_force_density"]:11.2e}')
    if failures:
        print(f'\n{len(failures)} case(s) skipped:')
        for lbl, err in failures:
            print(f'  {lbl}: {err}')

    if args.output and summaries:
        with open(args.output, 'w', newline='') as f:
            wtr = csv.DictWriter(f, fieldnames=list(summaries[0]))
            wtr.writeheader()
            wtr.writerows(summaries)
        print(f'\nsummary -> {args.output}')
    if args.per_snapshot and snaps:
        with open(args.per_snapshot, 'w', newline='') as f:
            wtr = csv.DictWriter(f, fieldnames=list(snaps[0]))
            wtr.writeheader()
            wtr.writerows(snaps)
        print(f'per-snapshot -> {args.per_snapshot}')


if __name__ == '__main__':
    main()
