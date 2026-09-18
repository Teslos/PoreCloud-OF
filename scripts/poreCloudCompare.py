#!/usr/bin/env python3
"""Compare many poreCloud runs on the same question poreCloudReport answers
for one: does the Leenov-Kolin exclusion force help push bubbles out of the
melt pool?

poreCloudReport.py already has every piece this needs - CSV loading, track
building, the buoyancy sign fix, and expulsion_budget(), which is the
decisive statistic.  This script does not reimplement any of that; it calls
poreCloudReport once per case, ranks the results, and adds a cross-case table
and PDF.  See poreCloudReport.py's module docstring for the unit contract
(CSV columns -> track dict) that expulsion_budget() relies on.

Usage examples:

    # explicit label=path pairs (path may be a CSV or a case dir containing
    # poreCloud_pores.csv)
    python3 poreCloudCompare.py \\
        control=/home/ivt/results/poreCloud-nofield/poreCloud_pores.csv \\
        dc-x=/home/ivt/results/poreCloud-dc-x \\
        azimuthal-0.2T=/home/ivt/results/poreCloud-0.2T-azimuthal-18rank \\
        -o compare.pdf

    # scan a directory of case directories, one poreCloud_pores.csv each,
    # label = subdirectory name
    python3 poreCloudCompare.py --dir /home/ivt/results -o compare.pdf

    # mix: everything under --dir plus one extra explicit case
    python3 poreCloudCompare.py --dir /home/ivt/results extra=/tmp/foo.csv

A case is treated as the no-field control if its label contains (after
stripping spaces/hyphens/underscores and lower-casing) 'control', 'nofield'
or 'baseline'; pass --control LABEL to force a specific one.  A missing or
empty CSV is reported as such (both in the table and the PDF) rather than
raising - the rest of the comparison still runs.
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import poreCloudReport as pcr


# ── discovering cases ────────────────────────────────────────────────────────

def discover_cases(dir_arg, explicit_specs):
    """Return an ordered list of (label, csv_path).

    --dir contributes one entry per subdirectory (label = dir name, path =
    <subdir>/poreCloud_pores.csv, whether or not it exists - existence is
    checked later so a missing run still shows up in the comparison).
    Explicit CASE args are label=path or a bare path (label defaults to the
    path's stem/dirname); a directory path is expanded to its
    poreCloud_pores.csv.  Explicit entries are appended after --dir entries,
    in the order given, and override a --dir entry of the same label.
    """
    cases = {}

    if dir_arg:
        d = Path(dir_arg)
        if not d.is_dir():
            raise SystemExit(f'--dir {dir_arg}: not a directory')
        for sub in sorted(p for p in d.iterdir() if p.is_dir()):
            cases[sub.name] = sub / 'poreCloud_pores.csv'

    for spec in explicit_specs:
        if '=' in spec:
            label, path_s = spec.split('=', 1)
        else:
            path_s, label = spec, None
        p = Path(path_s)
        if p.is_dir():
            p = p / 'poreCloud_pores.csv'
        if label is None:
            label = p.parent.name if p.name == 'poreCloud_pores.csv' else p.stem
        cases[label] = p

    if not cases:
        raise SystemExit('no cases given: pass --dir and/or label=path arguments')
    return list(cases.items())


# ── per-case loading ─────────────────────────────────────────────────────────

def load_case(label, path, min_len=4, cell_volume=pcr.DEFAULT_CELL_VOLUME,
             coverage_floor=pcr.DEFAULT_COVERAGE_FLOOR):
    """Load one case and reduce it to the row this comparison ranks on.

    Never raises for a missing/empty/malformed CSV - that is reported through
    the 'status' field ('ok', 'missing', 'empty', 'error') so one bad case
    cannot take down the whole comparison.  Only 'ok' rows carry the numeric
    fields the table and plots use.
    """
    path = Path(path)
    row = {'label': label, 'path': str(path), 'status': 'ok', 'error': None}

    if not path.exists():
        row['status'] = 'missing'
        row['error'] = f'no such file: {path}'
        return row
    if path.stat().st_size == 0:
        row['status'] = 'empty'
        row['error'] = 'file is empty'
        return row

    try:
        rows = pcr.load_rows(path)
    except SystemExit as e:
        # load_rows raises SystemExit for "no data rows" and for the
        # pre-Coverage schema; either way there is nothing to compare.
        row['status'] = 'empty'
        row['error'] = str(e)
        return row
    except (OSError, csv.Error) as e:
        row['status'] = 'error'
        row['error'] = f'{type(e).__name__}: {e}'
        return row

    try:
        tracks = pcr.build_tracks(rows, cell_volume=cell_volume, min_len=min_len)
        if not tracks:
            row['status'] = 'empty'
            row['error'] = f'no track reached {min_len} samples'
            return row
        budget = pcr.expulsion_budget(tracks)
        cov = pcr.coverage_stats(rows, floor=coverage_floor)
    except Exception as e:  # noqa: BLE001 - one bad case must not abort the run
        row['status'] = 'error'
        row['error'] = f'{type(e).__name__}: {e}'
        return row

    row.update({
        'budget': budget,
        'coverage': cov,
        'n_rows': len(rows),
        'n_bubbles': budget['n_bubbles'],
        'em_med': budget['em_median_mag'],
        'bu_med': budget['bu_median_mag'],
        'ratio': budget['ratio_median'],
        'up_frac': budget['em_up_fraction_bubbles'],
        'rose_frac': budget['rose_fraction'],
        'median_rise': budget['median_rise_um'],
        'reached_surface': budget['reached_surface'],
    })
    return row


def load_all(cases, **kw):
    return [load_case(label, path, **kw) for label, path in cases]


# ── ranking and baseline selection ──────────────────────────────────────────

def _rank_key(r):
    """Most-promising-first: does the exclusion force lift bubbles out?

    Primary key is the direction statistic the science question is actually
    decided on - the fraction of bubbles with a net upward EM impulse
    (poreCloudReport's own verdict threshold, 0.6/0.4, is built on exactly
    this number). Ties break on the outcome (did they reach the surface, did
    they rise) and finally on how large the EM force is relative to
    buoyancy - magnitude alone says nothing about whether it helps, so it is
    the last tiebreaker, not the first key.
    """
    return (-r['up_frac'], -r['reached_surface'], -r['rose_frac'], -r['ratio'])


def rank_cases(results):
    """Split into (ranked 'ok' cases, other-status cases, input order)."""
    valid = [r for r in results if r['status'] == 'ok']
    invalid = [r for r in results if r['status'] != 'ok']
    return sorted(valid, key=_rank_key), invalid


def _norm_label(s):
    return s.lower().replace('-', '').replace('_', '').replace(' ', '')


CONTROL_KEYWORDS = ('control', 'nofield', 'baseline')


def identify_control(results, explicit=None):
    """Pick the no-field control case, or None if none can be identified.

    explicit, when given, must match a case label exactly (raises otherwise -
    a typo silently falling back to auto-detection would be worse than
    failing loudly). Auto-detection matches CONTROL_KEYWORDS against the
    label with separators stripped, so 'no-field', 'No Field', 'nofield' and
    'NOFIELD-control' all match; on more than one match the first (input
    order) is used and a warning is printed - silently picking one without
    saying so is exactly the kind of thing that produces a plausible-looking
    but wrong table.
    """
    if explicit is not None:
        for r in results:
            if r['label'] == explicit:
                return r
        raise SystemExit(f'--control {explicit!r} matches no case label '
                          f'(have: {", ".join(r["label"] for r in results)})')

    matches = [r for r in results if r['status'] == 'ok'
              and any(k in _norm_label(r['label']) for k in CONTROL_KEYWORDS)]
    if len(matches) > 1:
        print(f"warning: multiple candidate control cases "
              f"({', '.join(m['label'] for m in matches)}); using "
              f"{matches[0]['label']!r} - pass --control to disambiguate",
              file=sys.stderr)
    return matches[0] if matches else None


# ── table ────────────────────────────────────────────────────────────────────

def format_table(ranked, invalid, control):
    """Render the ranked comparison table as a list of text lines."""
    lines = []
    W = 130
    lines.append('=' * W)
    lines.append('POREBLOUD CROSS-CASE COMPARISON - ranked most-promising-for-'
                 'expulsion first')
    lines.append('=' * W)
    if control is not None:
        lines.append(f"control (no-field baseline): {control['label']!r}")
    else:
        lines.append('control: none identified (no label matched control/'
                     'nofield/baseline; pass --control to set one)')
    lines.append('')

    hdr = (f"{'#':>2} {'label':<22} {'n_bub':>6} {'med|EM|(N)':>11} "
          f"{'med|buoy|(N)':>12} {'ratio':>7} {'%up':>6} {'%rose':>6} "
          f"{'med rise(um)':>13} {'%surf':>6} {'d%up vs ctrl':>13}")
    lines.append(hdr)
    lines.append('-' * len(hdr))

    for i, r in enumerate(ranked, 1):
        tag = ' *' if control is not None and r is control else ''
        if control is not None and r is not control:
            dup = f"{(r['up_frac'] - control['up_frac']) * 100:+.1f}pp"
        else:
            dup = '-' if control is None else '(baseline)'
        lines.append(
            f"{i:>2} {r['label']+tag:<22} {r['n_bubbles']:>6} "
            f"{r['em_med']:>11.3e} {r['bu_med']:>12.3e} {r['ratio']:>6.1f}x "
            f"{r['up_frac']*100:>5.1f}% {r['rose_frac']*100:>5.1f}% "
            f"{r['median_rise']:>+12.1f} {r['reached_surface']*100:>5.1f}% "
            f"{dup:>13}"
        )

    if invalid:
        lines.append('')
        lines.append('unavailable (excluded from ranking):')
        for r in invalid:
            lines.append(f"   {r['label']:<22} [{r['status']}] {r['error']}")

    lines.append('=' * W)
    return lines


def print_table(ranked, invalid, control):
    for line in format_table(ranked, invalid, control):
        print(line)


# ── PDF ──────────────────────────────────────────────────────────────────────

def make_comparison_pdf(ranked, invalid, control, colors, output_path):
    """Multi-page PDF: the table, then three cross-case plots."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    mpl.rcParams.update({
        'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11,
        'legend.fontsize': 9, 'figure.facecolor': 'white',
        'axes.facecolor': '#f8f8f8', 'axes.grid': True,
        'grid.alpha': 0.3, 'grid.linewidth': 0.5,
    })

    labels = [r['label'] for r in ranked]
    n = len(ranked)
    bar_colors = [colors[i % len(colors)] for i in range(n)]
    control_i = next((i for i, r in enumerate(ranked) if r is control), None)

    def finish(fig):
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    with PdfPages(output_path) as pdf:

        # Page 0: the table, as text - keeps the PDF self-contained.
        fig = plt.figure(figsize=(13, max(4, 0.35 * (n + len(invalid) + 8))))
        fig.text(0.02, 0.98, '\n'.join(format_table(ranked, invalid, control)),
                 family='monospace', fontsize=7.5, va='top', ha='left')
        plt.axis('off')
        pdf.savefig(fig)
        plt.close(fig)

        if n == 0:
            finish_fig = plt.figure(figsize=(8, 3))
            plt.text(0.5, 0.5, 'no usable cases to plot', ha='center', va='center')
            plt.axis('off')
            pdf.savefig(finish_fig)
            plt.close(finish_fig)
            return

        # Page 1: % of bubbles pushed upward, per case, ranked order.
        fig, ax = plt.subplots(figsize=(max(7, 0.9 * n), 5))
        up = [r['up_frac'] * 100 for r in ranked]
        bars = ax.bar(labels, up, color=bar_colors)
        if control_i is not None:
            bars[control_i].set_edgecolor('black')
            bars[control_i].set_linewidth(2.0)
            ax.axhline(up[control_i], color='k', ls='--', lw=1,
                      label=f"control ({control['label']})")
            ax.legend()
        ax.axhline(50, color='gray', lw=0.8, ls=':')
        ax.set_ylabel('bubbles with net UPWARD exclusion impulse  (%)')
        ax.set_title('Does the exclusion force push bubbles toward the surface?')
        plt.setp(ax.get_xticklabels(), rotation=35, ha='right')
        finish(fig)

        # Page 2: force-magnitude ratio, per case.
        fig, ax = plt.subplots(figsize=(max(7, 0.9 * n), 5))
        ratios = [r['ratio'] for r in ranked]
        bars = ax.bar(labels, ratios, color=bar_colors)
        if control_i is not None:
            bars[control_i].set_edgecolor('black')
            bars[control_i].set_linewidth(2.0)
        ax.axhline(1.0, color='gray', lw=0.8, ls=':', label='|EM| = |buoyancy|')
        ax.set_yscale('log')
        ax.set_ylabel('median |F_exclusion| / |F_buoyancy|')
        ax.set_title('How large is the exclusion force relative to buoyancy?')
        ax.legend()
        plt.setp(ax.get_xticklabels(), rotation=35, ha='right')
        finish(fig)

        # Page 3: net vertical displacement distributions, small multiples.
        ncols = min(3, n)
        nrows = -(-n // ncols)
        all_rise = np.concatenate([r['budget']['bubble']['rise'] for r in ranked])
        finite = all_rise[np.isfinite(all_rise)]
        lim = np.percentile(np.abs(finite), 98) if len(finite) else 1.0
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(4.2 * ncols, 3.2 * nrows), squeeze=False)
        fig.suptitle('Net vertical displacement per bubble '
                     '(+ = toward the free surface)', fontsize=13)
        for i, r in enumerate(ranked):
            ax = axes[i // ncols][i % ncols]
            rise = r['budget']['bubble']['rise']
            ax.hist(np.clip(rise, -lim, lim), bins=40, color=bar_colors[i], alpha=.85)
            ax.axvline(0, color='k', lw=1)
            title = r['label'] + (' (control)' if r is control else '')
            ax.set_title(f"{title}\n{r['rose_frac']*100:.0f}% rose, "
                        f"median {r['median_rise']:+.0f} um", fontsize=9)
        for j in range(n, nrows * ncols):
            axes[j // ncols][j % ncols].axis('off')
        finish(fig)

        # Page 4: overview - ratio vs direction, one point per case.
        fig, ax = plt.subplots(figsize=(7.5, 6))
        for i, r in enumerate(ranked):
            ax.scatter(r['ratio'], r['up_frac'] * 100, s=90, color=bar_colors[i],
                      edgecolor='black' if r is control else 'none',
                      linewidth=1.6, zorder=3)
            ax.annotate(r['label'], (r['ratio'], r['up_frac'] * 100),
                       xytext=(5, 5), textcoords='offset points', fontsize=8)
        ax.axhline(50, color='gray', lw=0.8, ls=':')
        ax.axvline(1, color='gray', lw=0.8, ls=':')
        ax.set_xscale('log')
        ax.set_xlabel('median |F_exclusion| / |F_buoyancy|')
        ax.set_ylabel('bubbles with net UPWARD exclusion impulse  (%)')
        ax.set_title('Magnitude vs direction - large does not mean helpful')
        finish(fig)


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cases', nargs='*',
                    help="label=path pairs (path = a poreCloud_pores.csv or a "
                         "case dir containing one); a bare path is also "
                         "accepted and gets an auto label")
    ap.add_argument('--dir', help='directory of case subdirectories to scan '
                                 '(label = subdir name)')
    ap.add_argument('-o', '--output', type=Path, default=Path('poreCloudCompare.pdf'),
                    help='PDF path (default: poreCloudCompare.pdf)')
    ap.add_argument('--control', help='label of the no-field control case '
                                      '(auto-detected if omitted)')
    ap.add_argument('--report-lib', help='directory holding analyse_pore_forces.py')
    ap.add_argument('--min-len', type=int, default=4, help='min samples per track')
    ap.add_argument('--cell-volume', type=float, default=pcr.DEFAULT_CELL_VOLUME)
    ap.add_argument('--coverage-floor', type=float, default=pcr.DEFAULT_COVERAGE_FLOOR)
    args = ap.parse_args()

    cases = discover_cases(args.dir, args.cases)
    results = load_all(cases, min_len=args.min_len, cell_volume=args.cell_volume,
                       coverage_floor=args.coverage_floor)

    ranked, invalid = rank_cases(results)
    control = identify_control(results, args.control)
    print_table(ranked, invalid, control)

    apf = pcr._import_report_layer(args.report_lib)
    make_comparison_pdf(ranked, invalid, control, apf.PUB_COLORS, args.output)
    print(f'\ncomparison -> {args.output}  '
         f'({len(ranked)} usable of {len(results)} cases)')


if __name__ == '__main__':
    main()
