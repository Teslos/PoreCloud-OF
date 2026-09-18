#!/usr/bin/env python3
"""Tests for poreCloudFieldSurvey.

Run:  python3 scripts/test_poreCloudFieldSurvey.py
      (bare runner, no pytest dependency - matches test_poreCloudReport.py)

This tool exists because trajectory replay on finished runs does not work
(see docs/specs.md); it answers the same "does the force point up?" question
from a single instantaneous field instead. Two things would make it silently
wrong rather than loudly wrong, and both get their own test:

  - the ASCII/BINARY auto-detection in read_field() could return the wrong
    numbers without raising (see docs on the finiteness/magnitude guard);
  - the sign flip from LorentzForce to the exclusion force could be backwards,
    which would make the tool report the exact opposite physics.

Everything else (masking, small-pool guard, zero-field NaN convention,
force-weighting) gets one hand-checkable case each.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import poreCloudFieldSurvey as pcfs


# ── OpenFOAM field-file fixtures ─────────────────────────────────────────────
#
# read_field() only cares about the internalField block, so these fixtures
# skip the rest of a real FoamFile header - a minimal preamble is included
# only for realism, not because the parser needs it.

def _ascii_scalar_text(values):
    body = '\n'.join(f'{v:.12g}' for v in values)
    return (
        'FoamFile\n{\n    format ascii;\n    class volScalarField;\n}\n\n'
        f'internalField   nonuniform List<scalar> \n{len(values)}\n(\n{body}\n)\n;\n\n'
        'boundaryField\n{\n}\n'
    )


def _ascii_vector_text(vectors):
    body = '\n'.join(f'({v[0]:.12g} {v[1]:.12g} {v[2]:.12g})' for v in vectors)
    return (
        'FoamFile\n{\n    format ascii;\n    class volVectorField;\n}\n\n'
        f'internalField   nonuniform List<vector> \n{len(vectors)}\n(\n{body}\n)\n;\n\n'
        'boundaryField\n{\n}\n'
    )


def _uniform_scalar_text(value):
    return f'FoamFile\n{{\n}}\ninternalField   uniform {value};\n'


def _uniform_vector_text(vec):
    return f'FoamFile\n{{\n}}\ninternalField   uniform ({vec[0]} {vec[1]} {vec[2]});\n'


def _write_binary_scalar(path, values):
    """Real little-endian float64 bytes, no newline between the '(' and the data."""
    n = len(values)
    header = (
        'FoamFile\n{\n    format binary;\n    class volScalarField;\n}\n\n'
        f'internalField   nonuniform List<scalar>\n{n}\n('
    ).encode('ascii')
    data = np.asarray(values, dtype='<f8').tobytes()
    with open(path, 'wb') as f:
        f.write(header + data + b')\n;\n\nboundaryField\n{\n}\n')


def _write_binary_vector(path, vectors):
    n = len(vectors)
    header = (
        'FoamFile\n{\n    format binary;\n    class volVectorField;\n}\n\n'
        f'internalField   nonuniform List<vector>\n{n}\n('
    ).encode('ascii')
    data = np.asarray(vectors, dtype='<f8').tobytes()
    with open(path, 'wb') as f:
        f.write(header + data + b')\n;\n\nboundaryField\n{\n}\n')


def _write(path, text):
    with open(path, 'w') as f:
        f.write(text)


def _make_snapshot(tmp_path, t='0', alpha=None, eps=None, lf=None):
    """Write a case/<t>/{alpha.metal,epsilon1,LorentzForce} snapshot.

    alpha/eps/lf are ('ascii'|'binary'|'uniform', values) tuples; ascii/binary
    write nonuniform fields, uniform writes a single broadcast value.
    """
    d = tmp_path / t
    d.mkdir(parents=True, exist_ok=True)
    _write_scalar(d / 'alpha.metal', *alpha)
    _write_scalar(d / 'epsilon1', *eps)
    _write_vector(d / 'LorentzForce', *lf)
    return tmp_path, t


def _write_scalar(path, kind, values):
    if kind == 'ascii':
        _write(path, _ascii_scalar_text(values))
    elif kind == 'binary':
        _write_binary_scalar(path, values)
    elif kind == 'uniform':
        _write(path, _uniform_scalar_text(values))
    else:
        raise ValueError(kind)


def _write_vector(path, kind, values):
    if kind == 'ascii':
        _write(path, _ascii_vector_text(values))
    elif kind == 'binary':
        _write_binary_vector(path, values)
    elif kind == 'uniform':
        _write(path, _uniform_vector_text(values))
    else:
        raise ValueError(kind)


# ── read_field: ascii ────────────────────────────────────────────────────────

def test_ascii_nonuniform_scalar_parses():
    import tempfile
    p = Path(tempfile.mkdtemp()) / 'alpha.metal'
    vals = [0.1, 0.25, 0.333333, 1.5, -2.75, 0.0]
    _write(p, _ascii_scalar_text(vals))
    a = pcfs.read_field(p, vec=False)
    assert np.allclose(a, vals)


def test_ascii_nonuniform_vector_parses():
    import tempfile
    p = Path(tempfile.mkdtemp()) / 'LorentzForce'
    vecs = [[0, 0, 0], [1, 2, 3], [-1.5, -2.5, -3.5], [0.5, 0.25, -0.75]]
    _write(p, _ascii_vector_text(vecs))
    a = pcfs.read_field(p, vec=True)
    assert a.shape == (4, 3)
    assert np.allclose(a, vecs)


# ── read_field: binary ───────────────────────────────────────────────────────

def test_binary_nonuniform_scalar_parses():
    import tempfile
    p = Path(tempfile.mkdtemp()) / 'epsilon1'
    vals = [1.234567, -9870.0, 0.0, 42.0, -1e-6, 3.14159265358979]
    _write_binary_scalar(p, vals)
    a = pcfs.read_field(p, vec=False)
    assert np.allclose(a, vals, rtol=0, atol=0)     # exact float64 round-trip


def test_binary_nonuniform_vector_parses():
    import tempfile
    p = Path(tempfile.mkdtemp()) / 'LorentzForce'
    vecs = [[1.1, -2.2, 3.3], [0.0, 0.0, 0.0], [-123456.789, 5e-8, -1.0]]
    _write_binary_vector(p, vecs)
    a = pcfs.read_field(p, vec=True)
    assert a.shape == (3, 3)
    assert np.allclose(a, vecs, rtol=0, atol=0)


def test_ascii_data_is_not_misread_as_binary():
    """The finiteness/magnitude guard exists because interpreting an ASCII
    text block as raw float64 does not raise - it returns numbers that look
    like data. If the guard ever failed to fall back to the ascii parser,
    this would come back as some large/NaN/inf garbage instead of the exact
    values below, so an exact match is sufficient to prove the ascii path
    was actually taken.

    NOTE ON VALUE CHOICE: this guard is weaker than it looks (see the bug
    noted in the module docstring/report - not fixed here per instructions).
    Reinterpreting short ASCII digit runs as float64 often lands on a small,
    finite, in-range number purely because ASCII byte values (0x30-0x39)
    happen to produce small exponents - e.g. [0.123456789, 0.987654321, 0.5,
    1e-06, 0.999999] round-trips as ~1e-33..1e-259 "data" that still passes
    `isfinite & abs<1e30`, and read_field silently returns it. The values
    below are chosen (and verified below) to actually exercise a case the
    guard catches, matching what happens on this project's real ASCII field
    files, where an out-of-[0,1]-range or large-magnitude value is always
    present somewhere in the block."""
    import tempfile
    p = Path(tempfile.mkdtemp()) / 'alpha.metal'
    vals = [0.1, 0.25, 0.333333, 1.5, -2.75, 0.0]
    _write(p, _ascii_scalar_text(vals))
    a = pcfs.read_field(p, vec=False)
    assert a.shape == (6,)
    assert np.allclose(a, vals)
    # and, independently: the raw text bytes really do fail the guard here,
    # so the fallback branch is not merely untested dead code for this input.
    raw = open(p, 'rb').read()
    start = raw.index(b'(\n') + 1
    block = raw[start:start + 8 * len(vals)]
    reinterpreted = np.frombuffer(block, dtype='<f8')
    assert not (np.isfinite(reinterpreted).all() and np.abs(reinterpreted).max() < 1e30)


# ── read_field: uniform ──────────────────────────────────────────────────────

def test_uniform_scalar_parses():
    import tempfile
    p = Path(tempfile.mkdtemp()) / 'alpha.metal'
    _write(p, _uniform_scalar_text(0.75))
    a = pcfs.read_field(p, vec=False)
    assert np.allclose(a, [0.75])


def test_uniform_vector_parses():
    import tempfile
    p = Path(tempfile.mkdtemp()) / 'LorentzForce'
    _write(p, _uniform_vector_text((1.0, -2.0, 3.0)))
    a = pcfs.read_field(p, vec=True)
    assert a.shape == (1, 3)
    assert np.allclose(a, [[1.0, -2.0, 3.0]])


def test_uniform_fields_broadcast_in_snapshot():
    """A uniform alpha/epsilon1 (single value) must broadcast to match a
    nonuniform LorentzForce of n cells, not truncate the snapshot to n=1."""
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    n = 150
    lf = [[0.0, -10.0, 0.0]] * n
    _make_snapshot(tmp, alpha=('uniform', 1.0), eps=('uniform', 1.0), lf=('ascii', lf))
    r = pcfs.survey_snapshot(str(tmp), '0')
    assert r is not None
    assert r['n_liquid'] == n


# ── liquid mask ───────────────────────────────────────────────────────────────

def test_liquid_mask_requires_both_alpha_and_epsilon_above_half():
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    # 120 cells: both fields liquid: alpha>0.5 and eps>0.5
    #  80 cells: alpha liquid, eps not
    #  80 cells: eps liquid, alpha not
    # 120 cells: neither
    # (the "both" group is kept >=100 so this exercises the mask, not the
    # separate small-pool guard tested below)
    alpha = [0.9] * 120 + [0.9] * 80 + [0.1] * 80 + [0.1] * 120
    eps = [0.9] * 120 + [0.1] * 80 + [0.9] * 80 + [0.1] * 120
    n = len(alpha)
    lf = [[0.0, -1.0, 0.0]] * n
    _make_snapshot(tmp, alpha=('ascii', alpha), eps=('ascii', eps), lf=('ascii', lf))
    r = pcfs.survey_snapshot(str(tmp), '0')
    assert r is not None
    assert r['n_liquid'] == 120


# ── the sign convention ──────────────────────────────────────────────────────

def test_lorentz_force_down_is_reported_as_exclusion_up():
    """This is the one that would silently invert every conclusion in the
    library if it were wrong. F_exclusion = -1.5 * V * LorentzForce, so a
    Lorentz body force pointing in -y (into the pool) must be reported as
    the exclusion force pointing +y (toward the free surface, frac_up ~ 1),
    not the other way around."""
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    n = 150
    alpha = [1.0] * n
    eps = [1.0] * n
    lf = [[0.0, -1000.0, 0.0]] * n     # Lorentz force points straight down (-y)
    _make_snapshot(tmp, alpha=('ascii', alpha), eps=('ascii', eps), lf=('ascii', lf))
    r = pcfs.survey_snapshot(str(tmp), '0')
    assert r is not None
    assert r['frac_up'] == 1.0, 'LorentzForce in -y must yield an exclusion force reported as fully up'
    assert r['median_ef_y'] > 0


def test_lorentz_force_up_is_reported_as_exclusion_down():
    """The mirror case, to pin the convention from both sides."""
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    n = 150
    alpha = [1.0] * n
    eps = [1.0] * n
    lf = [[0.0, 1000.0, 0.0]] * n       # Lorentz force points straight up (+y)
    _make_snapshot(tmp, alpha=('ascii', alpha), eps=('ascii', eps), lf=('ascii', lf))
    r = pcfs.survey_snapshot(str(tmp), '0')
    assert r is not None
    assert r['frac_up'] == 0.0, 'LorentzForce in +y must yield an exclusion force reported as fully down'
    assert r['median_ef_y'] < 0


# ── zero-field control ───────────────────────────────────────────────────────

def test_zero_force_field_returns_nan_not_zero():
    """No field -> no direction to ask about. frac_up must be NaN so the
    control is not mistaken for the worst-performing (0%) configuration."""
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    n = 150
    alpha = [1.0] * n
    eps = [1.0] * n
    lf = [[0.0, 0.0, 0.0]] * n
    _make_snapshot(tmp, alpha=('ascii', alpha), eps=('ascii', eps), lf=('ascii', lf))
    r = pcfs.survey_snapshot(str(tmp), '0')
    assert r is not None
    assert r['n_liquid'] == n
    assert math.isnan(r['frac_up'])
    assert math.isnan(r['frac_up_force_weighted'])


# ── small-pool guard ──────────────────────────────────────────────────────────

def test_small_pool_returns_none():
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    n = 200
    n_liquid = 50                          # fewer than the 100-cell floor
    alpha = [1.0] * n_liquid + [0.0] * (n - n_liquid)
    eps = [1.0] * n_liquid + [0.0] * (n - n_liquid)
    lf = [[0.0, -1.0, 0.0]] * n
    _make_snapshot(tmp, alpha=('ascii', alpha), eps=('ascii', eps), lf=('ascii', lf))
    r = pcfs.survey_snapshot(str(tmp), '0')
    assert r is None


def test_pool_just_at_the_floor_is_kept():
    """100 liquid cells is the floor, not the cutoff below it."""
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    n = 200
    n_liquid = 100
    alpha = [1.0] * n_liquid + [0.0] * (n - n_liquid)
    eps = [1.0] * n_liquid + [0.0] * (n - n_liquid)
    lf = [[0.0, -1.0, 0.0]] * n
    _make_snapshot(tmp, alpha=('ascii', alpha), eps=('ascii', eps), lf=('ascii', lf))
    r = pcfs.survey_snapshot(str(tmp), '0')
    assert r is not None
    assert r['n_liquid'] == 100


# ── force-weighting ───────────────────────────────────────────────────────────

def test_force_weighted_fraction_differs_from_plain_when_few_cells_dominate():
    """190 liquid cells carry a small downward force (0.01 N/m3, "down");
    10 carry a huge upward force (1000 N/m3, "up"). Plain frac_up counts
    cells (10/200 = 5%); force-weighted counts force (the 10 upward cells
    carry ~99.98% of the total |F|), so the two must disagree sharply."""
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    n_down, n_up = 190, 10
    n = n_down + n_up
    alpha = [1.0] * n
    eps = [1.0] * n
    lf = [[0.0, 0.01, 0.0]] * n_down + [[0.0, -1000.0, 0.0]] * n_up
    _make_snapshot(tmp, alpha=('ascii', alpha), eps=('ascii', eps), lf=('ascii', lf))
    r = pcfs.survey_snapshot(str(tmp), '0')
    assert r is not None
    assert math.isclose(r['frac_up'], n_up / n, rel_tol=1e-9)
    expected_fw = (n_up * 1000.0) / (n_up * 1000.0 + n_down * 0.01)
    assert math.isclose(r['frac_up_force_weighted'], expected_fw, rel_tol=1e-6)
    assert r['frac_up_force_weighted'] > 0.99
    assert r['frac_up'] < 0.10
    assert r['frac_up_force_weighted'] - r['frac_up'] > 0.5


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
