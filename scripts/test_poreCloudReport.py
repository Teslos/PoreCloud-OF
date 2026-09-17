#!/usr/bin/env python3
"""Tests for poreCloudReport.

Run:  python3 -m pytest scripts/test_poreCloudReport.py -q
      (or plain `python3 scripts/test_poreCloudReport.py` for a bare runner)

The conversions are the whole point of the adapter, and every one of them is a
place where a wrong answer looks entirely plausible on a plot: metres read as
micrometres, a force read as a force density, a slip velocity with the sign
flipped.  Each gets a test with a hand-checkable number.
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import poreCloudReport as pcr

COLUMNS = [
    'Time', 'PoreID', 'IsKeyhole', 'Volume_m3', 'Cx_m', 'Cy_m', 'Cz_m',
    'Fx_N', 'Fy_N', 'Fz_N', 'BirthTime', 'BirthCx', 'BirthCy', 'BirthCz',
    'NCells', 'Ux', 'Uy', 'Uz', 'd_m', 'dOverDx', 'Active',
    'Uslipx', 'Uslipy', 'Uslipz', 'Fbuoyx', 'Fbuoyy', 'Fbuoyz',
    'T_K', 'epsilon1', 'Coverage',
]

DEFAULTS = {c: '0' for c in COLUMNS}
DEFAULTS.update({'Volume_m3': '1e-15', 'Active': '1', 'Coverage': '1',
                 'T_K': '1800', 'IsKeyhole': '0', 'NCells': '50',
                 'd_m': '3.5e-5', 'dOverDx': '4.2', 'epsilon1': '1'})


def make_rows(specs):
    """specs: list of dicts overriding DEFAULTS."""
    return [dict(DEFAULTS, **s) for s in specs]


def write_csv(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)


# ── unit conversions ─────────────────────────────────────────────────────────

def test_positions_convert_m_to_um():
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Cx_m': '1e-4',
         'Cy_m': '2e-4', 'Cz_m': '3e-4'} for i in range(4)
    ])
    td = pcr.build_tracks(rows)[0]
    assert np.allclose(td['x'], 100.0)      # 1e-4 m = 100 um
    assert np.allclose(td['depth'], 200.0)
    assert np.allclose(td['width'], 300.0)


def test_time_converts_s_to_us():
    rows = make_rows([{'PoreID': '1', 'Time': f'{i * 1e-6}'} for i in range(4)])
    td = pcr.build_tracks(rows)[0]
    assert np.allclose(td['t_us'], [0.0, 1.0, 2.0, 3.0])


def test_force_becomes_density_by_dividing_volume():
    """2e-9 N in a 1e-15 m3 bubble is 2e6 N/m3 - the report's axis unit."""
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Fx_N': '2e-9',
         'Volume_m3': '1e-15', 'Fbuoyy': '5e-10'} for i in range(4)
    ])
    td = pcr.build_tracks(rows)[0]
    assert np.allclose(td['EF'][:, 0], 2e6)
    assert np.allclose(td['LF'][:, 1], 5e5)


def test_carrier_velocity_subtracts_slip():
    """Uslip = U_parcel - Uc, so Uc = U - Uslip.  Sign errors land here."""
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Ux': '1.0', 'Uslipx': '0.25'}
        for i in range(4)
    ])
    td = pcr.build_tracks(rows)[0]
    assert np.allclose(td['U_pore'][:, 0], 0.75)


def test_size_is_volume_in_cells():
    cell_v = 2e-15
    rows = make_rows([{'PoreID': '1', 'Time': f'{i * 1e-6}',
                       'Volume_m3': '6e-15'} for i in range(4)])
    td = pcr.build_tracks(rows, cell_volume=cell_v)[0]
    assert np.allclose(td['size'], 3.0)


# ── kinematics ───────────────────────────────────────────────────────────────

def test_velocity_um_per_us_equals_m_per_s():
    """Moving 1 um per us is 1 m/s; the report's axes assume that identity."""
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Cx_m': f'{i * 1e-6}'}
        for i in range(5)
    ])
    td = pcr.build_tracks(rows)[0]
    assert np.allclose(td['vx'], 1.0)
    assert np.allclose(td['vmag'], 1.0)


def test_displacement_and_lifetime():
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Cx_m': f'{i * 1e-6}'}
        for i in range(5)
    ])
    td = pcr.build_tracks(rows)[0]
    assert np.allclose(td['displacement'], [4.0, 0.0, 0.0])
    assert math.isclose(td['lifetime_us'], 4.0)


def test_impulse_is_force_density_times_time():
    """Constant 1e6 N/m3 for 4 us -> 4 N*s/m3 over the 4 leading intervals."""
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Fx_N': '1e-9',
         'Volume_m3': '1e-15'} for i in range(5)
    ])
    td = pcr.build_tracks(rows)[0]
    assert math.isclose(td['EF_impulse'][0], 1e6 * 4e-6, rel_tol=1e-9)


# ── grouping and filtering ───────────────────────────────────────────────────

def test_tracks_group_by_poreid_not_position():
    """Persistent ids are why this adapter needs no nearest-neighbour matching."""
    rows = make_rows(
        [{'PoreID': '1', 'Time': f'{i * 1e-6}', 'Cx_m': '0'} for i in range(4)] +
        [{'PoreID': '2', 'Time': f'{i * 1e-6}', 'Cx_m': '0'} for i in range(4)]
    )
    assert len(pcr.build_tracks(rows)) == 2


def test_short_tracks_dropped():
    rows = make_rows([{'PoreID': '1', 'Time': f'{i * 1e-6}'} for i in range(3)])
    assert pcr.build_tracks(rows, min_len=4) == []
    assert len(pcr.build_tracks(rows, min_len=3)) == 1


def test_captured_pores_excluded_from_tracks():
    """A frozen pore holds position by design; its zero velocity is not physics."""
    rows = make_rows([{'PoreID': '1', 'Time': f'{i * 1e-6}', 'Active': '0'}
                      for i in range(4)])
    assert pcr.build_tracks(rows, active_only=True) == []
    assert len(pcr.build_tracks(rows, active_only=False)) == 1


def test_rows_out_of_time_order_are_sorted():
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Cx_m': f'{i * 1e-6}'}
        for i in [3, 0, 2, 1]
    ])
    td = pcr.build_tracks(rows)[0]
    assert np.all(np.diff(td['t_us']) > 0)
    assert np.allclose(td['x'], [0.0, 1.0, 2.0, 3.0])


def test_zero_volume_does_not_raise_or_poison_impulse():
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Fx_N': '1e-9',
         'Volume_m3': '0' if i == 1 else '1e-15'} for i in range(5)
    ])
    td = pcr.build_tracks(rows)[0]
    assert not np.isfinite(td['EF'][1, 0])          # that sample is NaN
    assert np.isfinite(td['EF_impulse'][0])         # but the sum survives


# ── summary statistics ───────────────────────────────────────────────────────

def test_coverage_stats_counts_below_floor():
    rows = make_rows([{'PoreID': str(i), 'Coverage': c}
                      for i, c in enumerate(['1.0', '0.95', '0.80', '0.40'])])
    cs = pcr.coverage_stats(rows, floor=0.89)
    assert cs['n'] == 4
    assert math.isclose(cs['frac_below_floor'], 0.5)     # 0.80 and 0.40
    assert math.isclose(cs['frac_below_half'], 0.25)     # 0.40
    assert math.isclose(cs['min'], 0.40)


def test_capture_stats_uses_final_snapshot_only():
    rows = make_rows(
        [{'PoreID': '1', 'Time': '1e-6', 'Active': '1'},
         {'PoreID': '2', 'Time': '1e-6', 'Active': '1'},
         {'PoreID': '1', 'Time': '2e-6', 'Active': '1'},
         {'PoreID': '2', 'Time': '2e-6', 'Active': '0'}]
    )
    ps = pcr.capture_stats(rows)
    assert ps['n_parcels'] == 2 and ps['n_active'] == 1 and ps['n_pores'] == 1


def test_force_discontinuity_separates_stencil_change_from_motion():
    """Construct the effect: same displacement, force steps only when coverage does."""
    rows = []
    for pid in range(120):
        for i in range(4):
            jumped = (pid % 2 == 0) and i >= 2
            rows.append(dict(
                DEFAULTS, PoreID=str(pid), Time=f'{i * 1e-6}',
                Cx_m=f'{i * 1e-7}',                       # same motion for all
                Coverage='0.5' if jumped else '1.0',      # half change stencil
                Fx_N='2e-9' if jumped else '1e-9',        # and those step 100%
            ))
    disc = pcr.force_discontinuity(pcr.build_tracks(rows), dx=1e-5)
    assert disc, 'expected at least one populated displacement bin'
    d = disc[0]
    assert d['median_steady'] < 1e-6        # unchanged stencil -> no step
    assert math.isclose(d['median_jump'], 1.0, rel_tol=1e-6)   # 100% step


# ── physical invariants ──────────────────────────────────────────────────────

def test_buoyancy_impulse_is_antiparallel_to_gravity():
    """Buoyancy must lie on the gravity axis alone.

    This is the invariant that caught a y/z transposition in the shared
    report layer's summary printer: a buoyancy impulse with a non-zero
    cross-gravity component is impossible, so it can only be a plumbing bug.
    """
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}',
         'Fbuoyx': '0', 'Fbuoyy': '-1.5e-9', 'Fbuoyz': '0'} for i in range(5)
    ])
    td = pcr.build_tracks(rows)[0]
    imp = td['LF_impulse']
    # The CSV here uses the pre-fix convention (buoyancy written along g);
    # build_tracks detects that and flips it, so the track must come back with
    # buoyancy opposing gravity - i.e. pointing +y, lifting the bubble.
    assert imp[1] > 0, f'buoyancy must oppose -y gravity, got {imp}'
    assert abs(imp[0]) < 1e-12 and abs(imp[2]) < 1e-12, \
        f'buoyancy has off-axis components: {imp}'


REAL_CSV = Path('/tmp/claude-1804035392/-home-ivt-mhd-laserbeamfoam-solver/'
                '8328c3eb-f374-424c-8b3d-b6018e3fd191/scratchpad/prod18/'
                'poreCloud_pores.csv')


def test_real_run_invariants():
    """Sanity-check the actual production CSV, when it is still on disk."""
    if not REAL_CSV.exists():
        print('    (skipped: production CSV not present)')
        return
    import itertools
    with open(REAL_CSV) as f:
        rows = list(itertools.islice(csv.DictReader(f), 20000))

    cov = np.array([float(r['Coverage']) for r in rows])
    assert cov.min() > 0 and cov.max() <= 1.0 + 1e-12, 'coverage out of (0,1]'

    vol = np.array([float(r['Volume_m3']) for r in rows])
    assert (vol > 0).all(), 'non-positive bubble volume'

    assert {int(float(r['Active'])) for r in rows} <= {0, 1}, 'Active not boolean'

    bx = np.array([float(r['Fbuoyx']) for r in rows])
    bz = np.array([float(r['Fbuoyz']) for r in rows])
    by = np.array([float(r['Fbuoyy']) for r in rows])
    assert (bx == 0).all() and (bz == 0).all(), 'buoyancy off the gravity axis'
    assert (by <= 0).all(), 'buoyancy not opposing -y gravity'

    t0 = min(float(r['Time']) for r in rows)
    first = [r['PoreID'] for r in rows if float(r['Time']) == t0]
    assert len(first) == len(set(first)), 'duplicate PoreID within one snapshot'


# ── buoyancy sign and expulsion budget ───────────────────────────────────────

def test_buoyancy_sign_detects_the_pre_fix_convention():
    """Pre-2026-09-17 CSVs wrote buoyancy along g; it must be flipped."""
    down = make_rows([{'PoreID': '1', 'Time': f'{i * 1e-6}', 'Fbuoyy': '-1e-9'}
                      for i in range(4)])
    up = make_rows([{'PoreID': '1', 'Time': f'{i * 1e-6}', 'Fbuoyy': '+1e-9'}
                    for i in range(4)])
    assert pcr.buoyancy_sign(down) == -1     # points along g -> needs flipping
    assert pcr.buoyancy_sign(up) == +1
    # and the flip must actually reach the track
    assert pcr.build_tracks(down)[0]['LF'][0, 1] > 0


def test_expulsion_budget_reads_a_helping_force_as_helping():
    """EM pushing +y on every bubble must come back as 100% upward."""
    rows = make_rows([
        {'PoreID': str(p), 'Time': f'{i * 1e-6}', 'Fy_N': '1e-9',
         'Fbuoyy': '1e-10', 'Cy_m': f'{(400 + i) * 1e-6}'}
        for p in range(5) for i in range(5)
    ])
    b = pcr.expulsion_budget(pcr.build_tracks(rows, buoy_sign=1))
    assert b['em_up_fraction_samples'] == 1.0
    assert b['em_up_fraction_bubbles'] == 1.0
    assert b['rose_fraction'] == 1.0
    assert b['median_rise_um'] > 0


def test_expulsion_budget_reads_an_opposing_force_as_opposing():
    rows = make_rows([
        {'PoreID': str(p), 'Time': f'{i * 1e-6}', 'Fy_N': '-1e-9',
         'Fbuoyy': '1e-10', 'Cy_m': f'{(400 - i) * 1e-6}'}
        for p in range(5) for i in range(5)
    ])
    b = pcr.expulsion_budget(pcr.build_tracks(rows, buoy_sign=1))
    assert b['em_up_fraction_bubbles'] == 0.0
    assert b['rose_fraction'] == 0.0
    assert b['impulse_ratio_median'] < 0      # opposes buoyancy


def test_expulsion_ratio_is_magnitude_not_signed():
    """A large downward EM force still counts as dominating the budget."""
    rows = make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Fy_N': '-1e-8', 'Fbuoyy': '1e-10'}
        for i in range(5)
    ])
    b = pcr.expulsion_budget(pcr.build_tracks(rows, buoy_sign=1))
    assert math.isclose(b['ratio_median'], 100.0, rel_tol=1e-6)
    assert b['em_exceeds_buoyancy'] == 1.0
    assert b['em_up_fraction_samples'] == 0.0   # ... while pointing the wrong way


# ── end-to-end ───────────────────────────────────────────────────────────────

def test_load_rows_rejects_pre_coverage_schema(tmp_path=None):
    import tempfile
    d = Path(tempfile.mkdtemp())
    p = d / 'old.csv'
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[c for c in COLUMNS if c != 'Coverage'])
        w.writeheader()
        w.writerow({c: '0' for c in COLUMNS if c != 'Coverage'})
    try:
        pcr.load_rows(p)
    except SystemExit as e:
        assert 'Coverage' in str(e)
    else:
        raise AssertionError('expected SystemExit on pre-Coverage schema')


def test_roundtrip_through_a_real_csv_file():
    import tempfile
    d = Path(tempfile.mkdtemp())
    p = d / 'poreCloud_pores.csv'
    write_csv(p, make_rows([
        {'PoreID': '1', 'Time': f'{i * 1e-6}', 'Cx_m': f'{i * 1e-6}'}
        for i in range(5)
    ]))
    tracks = pcr.build_tracks(pcr.load_rows(p))
    assert len(tracks) == 1 and np.allclose(tracks[0]['vx'], 1.0)


if __name__ == '__main__':
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith('test_') and callable(f)]
    failed = 0
    for name, fn in fns:
        try:
            fn()
            print(f'  PASS  {name}')
        except Exception as e:
            failed += 1
            print(f'  FAIL  {name}: {type(e).__name__}: {e}')
    print(f'\n{len(fns) - failed}/{len(fns)} passed')
    sys.exit(1 if failed else 0)
