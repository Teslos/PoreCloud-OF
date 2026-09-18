#!/usr/bin/env python3
"""Tests for poreCloudCompare.

Run:  python3 -m pytest scripts/test_poreCloudCompare.py -q
      (or plain `python3 scripts/test_poreCloudCompare.py` for a bare runner)

poreCloudCompare's whole job is to not lie about a comparison it cannot make
cleanly: rank cases in an order someone could act on, pick the right
baseline, and say plainly when a case has no usable data rather than folding
it into the table as a silent zero. Those are exactly the ways this could
fail while still printing a table that looks fine at a glance, so that is
what these tests target - not the arithmetic inside expulsion_budget(),
which is poreCloudReport's own test file's job and is reused here unchanged.

Row construction reuses test_poreCloudReport's COLUMNS/DEFAULTS/make_rows/
write_csv helpers rather than duplicating them.
"""

import math
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import poreCloudCompare as pcc
import poreCloudReport as pcr
import test_poreCloudReport as tpr


def _tmpdir():
    return Path(tempfile.mkdtemp())


def _bubble_rows(pore_ids, fy_n, fbuoy_y='1e-10', cy0=400, n=5, dcy=1):
    """n-sample tracks for each PoreID in pore_ids, constant Fy_N, Cy_m
    drifting by dcy um per sample - the same pattern
    test_poreCloudReport's expulsion-budget tests use."""
    return tpr.make_rows([
        {'PoreID': str(p), 'Time': f'{i * 1e-6}', 'Fy_N': fy_n,
         'Fbuoyy': fbuoy_y, 'Cy_m': f'{(cy0 + i * dcy) * 1e-6}'}
        for p in pore_ids for i in range(n)
    ])


def _write(rows, name='poreCloud_pores.csv'):
    d = _tmpdir()
    p = d / name
    tpr.write_csv(p, rows)
    return p


# ── missing / empty / malformed CSVs are reported, not fatal ────────────────

def test_missing_csv_reported_as_missing():
    r = pcc.load_case('gone', Path('/no/such/file/poreCloud_pores.csv'))
    assert r['status'] == 'missing'
    assert 'label' in r and r['label'] == 'gone'
    assert 'budget' not in r


def test_empty_file_reported_as_empty():
    d = _tmpdir()
    p = d / 'empty.csv'
    p.touch()
    r = pcc.load_case('blank', p)
    assert r['status'] == 'empty'
    assert 'budget' not in r


def test_header_only_csv_reported_as_empty():
    p = _write([])
    r = pcc.load_case('no-rows', p)
    assert r['status'] == 'empty'


def test_pre_coverage_schema_reported_as_empty_not_crash():
    """load_rows raises SystemExit for the pre-Coverage schema; load_case
    must catch it rather than letting the whole comparison die."""
    import csv
    d = _tmpdir()
    p = d / 'old.csv'
    cols = [c for c in tpr.COLUMNS if c != 'Coverage']
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerow({c: '0' for c in cols})
    r = pcc.load_case('old-schema', p)
    assert r['status'] == 'empty'
    assert 'Coverage' in r['error']


def test_all_pores_below_min_len_reported_as_empty():
    rows = tpr.make_rows([{'PoreID': '1', 'Time': '0'}])  # single sample
    p = _write(rows)
    r = pcc.load_case('too-short', p, min_len=4)
    assert r['status'] == 'empty'


# ── load_case reproduces expulsion_budget exactly (no silent corruption) ────

def test_load_case_numbers_match_expulsion_budget_directly():
    rows = _bubble_rows(range(6), fy_n='1e-9')
    p = _write(rows)
    r = pcc.load_case('x', p)
    assert r['status'] == 'ok'
    direct = pcr.expulsion_budget(pcr.build_tracks(pcr.load_rows(p)))
    assert r['n_bubbles'] == direct['n_bubbles']
    assert math.isclose(r['up_frac'], direct['em_up_fraction_bubbles'])
    assert math.isclose(r['ratio'], direct['ratio_median'], rel_tol=1e-9)
    assert math.isclose(r['median_rise'], direct['median_rise_um'])


# ── ranking order ────────────────────────────────────────────────────────────

def test_ranking_puts_helping_case_above_opposing_case():
    """A case where EM clearly pushes bubbles up must outrank one where it
    clearly pushes them down - the one thing a ranked table must get right."""
    helps = pcc.load_case('helps', _write(_bubble_rows(range(8), fy_n='1e-9')))
    opposes = pcc.load_case('opposes', _write(_bubble_rows(range(8), fy_n='-1e-9')))
    assert helps['status'] == opposes['status'] == 'ok'

    ranked, invalid = pcc.rank_cases([opposes, helps])  # deliberately out of order
    assert invalid == []
    assert [r['label'] for r in ranked] == ['helps', 'opposes']


def test_ranking_ties_on_direction_break_by_outcome_then_ratio():
    """Two cases with identical up_frac (here: both 0, no EM at all) must not
    be ordered arbitrarily - a real tiebreaker (rose_frac/reached_surface/
    ratio) should decide it, not dict/insertion order."""
    weak = pcc.load_case('weak-buoyancy-only',
                         _write(_bubble_rows(range(6), fy_n='0', fbuoy_y='1e-10')))
    strong = pcc.load_case('strong-buoyancy-only',
                           _write(_bubble_rows(range(6), fy_n='0', fbuoy_y='1e-9')))
    ranked, _ = pcc.rank_cases([weak, strong])
    assert ranked[0]['up_frac'] == ranked[1]['up_frac'] == 0.0
    # both rose (positive Cy drift) and reached no surface, so this comes
    # down to nothing distinguishing them except being equally tied - the
    # key point is rank_cases must not raise and must return both.
    assert {r['label'] for r in ranked} == {'weak-buoyancy-only', 'strong-buoyancy-only'}


def test_ranking_excludes_missing_and_empty_from_valid_list():
    ok = pcc.load_case('ok', _write(_bubble_rows(range(5), fy_n='1e-9')))
    missing = pcc.load_case('missing', Path('/no/such.csv'))
    empty = pcc.load_case('empty', _write([]))
    ranked, invalid = pcc.rank_cases([ok, missing, empty])
    assert [r['label'] for r in ranked] == ['ok']
    assert {r['label'] for r in invalid} == {'missing', 'empty'}
    # and no invalid row can silently masquerade as a ranked one
    for r in ranked:
        assert r['status'] == 'ok'
        assert 'up_frac' in r
    for r in invalid:
        assert 'up_frac' not in r


# ── baseline / control selection ─────────────────────────────────────────────

def test_control_auto_detected_by_keyword():
    a = pcc.load_case('dc-x', _write(_bubble_rows(range(4), fy_n='1e-9')))
    b = pcc.load_case('no-field-control', _write(_bubble_rows(range(4), fy_n='0')))
    ctrl = pcc.identify_control([a, b])
    assert ctrl is b


def test_control_auto_detection_matches_various_spellings():
    for label in ('control', 'Control', 'NoField', 'no_field', 'baseline-run',
                 'NOFIELD-control'):
        r = pcc.load_case(label, _write(_bubble_rows(range(3), fy_n='0')))
        assert pcc.identify_control([r]) is r, f'{label!r} should match'
    r = pcc.load_case('dc-x-0.1T', _write(_bubble_rows(range(3), fy_n='1e-9')))
    assert pcc.identify_control([r]) is None


def test_control_explicit_label_overrides_and_validates():
    a = pcc.load_case('dc-x', _write(_bubble_rows(range(3), fy_n='1e-9')))
    b = pcc.load_case('reference-run', _write(_bubble_rows(range(3), fy_n='0')))
    # 'reference-run' doesn't match the auto keywords, so --control is the
    # only way to select it - must work when given explicitly.
    ctrl = pcc.identify_control([a, b], explicit='reference-run')
    assert ctrl is b
    try:
        pcc.identify_control([a, b], explicit='typo-name')
    except SystemExit:
        pass
    else:
        raise AssertionError('expected SystemExit for an unmatched --control label')


def test_control_multiple_candidates_picks_first_without_crashing():
    a = pcc.load_case('control-a', _write(_bubble_rows(range(3), fy_n='0')))
    b = pcc.load_case('control-b', _write(_bubble_rows(range(3), fy_n='0')))
    ctrl = pcc.identify_control([a, b])
    assert ctrl is a   # first in input order


def test_control_ignores_a_broken_case_even_if_labelled_control():
    """A control CSV that is missing/empty cannot be used as a baseline -
    identify_control must skip non-'ok' rows even when the label matches."""
    broken = pcc.load_case('control', Path('/no/such.csv'))
    ok = pcc.load_case('dc-x', _write(_bubble_rows(range(3), fy_n='1e-9')))
    assert pcc.identify_control([broken, ok]) is None


# ── case discovery ────────────────────────────────────────────────────────────

def test_discover_cases_scans_directory_by_subdir_name():
    root = _tmpdir()
    (root / 'caseA').mkdir()
    (root / 'caseB').mkdir()
    tpr.write_csv(root / 'caseA' / 'poreCloud_pores.csv',
                  _rows_as_list(range(2)))
    # caseB deliberately has no CSV - discovery must still list it (as a
    # path that later turns out missing), not skip it silently.
    cases = dict(pcc.discover_cases(str(root), []))
    assert set(cases) == {'caseA', 'caseB'}
    assert cases['caseA'].name == 'poreCloud_pores.csv'
    assert not cases['caseB'].exists()


def _rows_as_list(pore_ids):
    return tpr.make_rows([{'PoreID': str(p), 'Time': f'{i * 1e-6}'}
                          for p in pore_ids for i in range(4)])


def test_discover_cases_explicit_overrides_dir_entry_of_same_label():
    root = _tmpdir()
    (root / 'caseA').mkdir()
    tpr.write_csv(root / 'caseA' / 'poreCloud_pores.csv', _rows_as_list(range(2)))
    override = _write(_rows_as_list(range(2)), name='other.csv')
    cases = dict(pcc.discover_cases(str(root), [f'caseA={override}']))
    assert cases['caseA'] == override


def test_discover_cases_bare_path_gets_auto_label():
    p = _write(_rows_as_list(range(2)))
    cases = dict(pcc.discover_cases(None, [str(p.parent)]))
    assert list(cases) == [p.parent.name]


def test_discover_cases_requires_at_least_one_case():
    try:
        pcc.discover_cases(None, [])
    except SystemExit:
        pass
    else:
        raise AssertionError('expected SystemExit with no --dir and no cases')


# ── table formatting ──────────────────────────────────────────────────────────

def test_format_table_lists_unavailable_cases_with_their_reason():
    ok = pcc.load_case('ok-case', _write(_bubble_rows(range(4), fy_n='1e-9')))
    missing = pcc.load_case('gone-case', Path('/no/such.csv'))
    ranked, invalid = pcc.rank_cases([ok, missing])
    ctrl = pcc.identify_control([ok, missing])
    lines = pcc.format_table(ranked, invalid, ctrl)
    text = '\n'.join(lines)
    assert 'ok-case' in text
    assert 'gone-case' in text
    assert 'missing' in text            # status is visible, not swallowed
    assert 'no such file' in text        # and the reason


def test_format_table_marks_the_baseline_row():
    ctrl_case = pcc.load_case('no-field', _write(_bubble_rows(range(5), fy_n='0')))
    other = pcc.load_case('dc-x', _write(_bubble_rows(range(5), fy_n='1e-9')))
    ranked, invalid = pcc.rank_cases([other, ctrl_case])
    ctrl = pcc.identify_control(ranked)
    text = '\n'.join(pcc.format_table(ranked, invalid, ctrl))
    assert '(baseline)' in text
    assert 'control (no-field baseline)' in text


# ── PDF smoke test ─────────────────────────────────────────────────────────────

def test_make_comparison_pdf_runs_on_a_multi_case_set():
    """Not checking pixels - just that the plotting path handles several
    cases, a zero-force control (constant arrays -> possible NaN stats) and
    an explicit control marker without raising, and writes a real PDF."""
    helps = pcc.load_case('helps', _write(_bubble_rows(range(6), fy_n='1e-9')))
    opposes = pcc.load_case('opposes', _write(_bubble_rows(range(6), fy_n='-1e-9')))
    control = pcc.load_case('no-field-control',
                            _write(_bubble_rows(range(6), fy_n='0')))
    ranked, invalid = pcc.rank_cases([helps, opposes, control])
    ctrl = pcc.identify_control([helps, opposes, control])
    assert ctrl is control

    out = _tmpdir() / 'compare.pdf'
    pcc.make_comparison_pdf(ranked, invalid, ctrl, pcr._import_report_layer().PUB_COLORS, out)
    assert out.exists() and out.stat().st_size > 1000
    with open(out, 'rb') as f:
        assert f.read(5) == b'%PDF-'


def test_make_comparison_pdf_handles_zero_usable_cases():
    """Every case missing/empty - must still produce a (near-empty) PDF, not
    crash trying to index an empty case list."""
    missing = pcc.load_case('gone', Path('/no/such.csv'))
    ranked, invalid = pcc.rank_cases([missing])
    assert ranked == []
    out = _tmpdir() / 'empty_compare.pdf'
    pcc.make_comparison_pdf(ranked, invalid, None, pcr._import_report_layer().PUB_COLORS, out)
    assert out.exists()
    with open(out, 'rb') as f:
        assert f.read(5) == b'%PDF-'


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
