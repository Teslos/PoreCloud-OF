#!/usr/bin/env python3
"""Tests for poreSeedingFromData.py.

The anchor test reproduces the weights currently in poreCloudProperties from
the PoreTracker summaries in data/. If the script and the dictionary ever
disagree, one of them has drifted and the seeding no longer means what its
comment says.
"""

import glob
import io
import os
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SCRIPT = os.path.join(HERE, "poreSeedingFromData.py")
sys.path.insert(0, HERE)

import poreSeedingFromData as M  # noqa: E402

# The six numbers presently in constant/poreCloudProperties.
IN_USE = [0.31264, 0.12417, 0.10200, 0.22395, 0.19290, 0.04435]
SUMMARIES = os.path.join(REPO, "data", "poretracker-summaries")


def run(args):
    buf = io.StringIO()
    with redirect_stdout(buf):
        M.main(args)
    return buf.getvalue()


def write_csv(rows, header):
    fh = tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, newline="")
    fh.write(",".join(header) + "\n")
    for r in rows:
        fh.write(",".join(str(x) for x in r) + "\n")
    fh.close()
    return fh.name


class PooledVoFWeights(unittest.TestCase):
    """The anchor: regenerate the in-use weights from the shipped data."""

    @classmethod
    def setUpClass(cls):
        files = sorted(glob.glob(os.path.join(SUMMARIES, "*_pores_summary.csv")))
        if not files:
            raise unittest.SkipTest(f"no summaries under {SUMMARIES}")
        # Pool the seven cases, keeping only non-keyhole pores, exactly as the
        # in-use weights were derived.
        import csv
        rows = []
        for f in files:
            with open(f, newline="") as fh:
                for r in csv.DictReader(fh):
                    if int(r["IsKeyhole"]) == 0:
                        rows.append((r["BirthCx_m"], r["BirthCz_m"], r["BirthTime_s"]))
        cls.path = write_csv(rows, ["x", "z", "t"])
        cls.n = len(rows)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.path)

    def test_pool_size(self):
        self.assertEqual(self.n, 452)

    def test_reproduces_in_use_weights(self):
        out = run([self.path, "--col", "x=x", "--col", "z=z", "--col", "t=t",
                   "--i-have-birth-data"])
        line = next(l for l in out.splitlines() if l.strip().startswith("weights"))
        got = [float(x) for x in line.split("(")[1].split(")")[0].split()]
        np.testing.assert_allclose(got, IN_USE, atol=5e-6)

    def test_one_pore_dropped_beyond_outer_edge(self):
        out = run([self.path, "--col", "x=x", "--col", "z=z", "--col", "t=t",
                   "--i-have-birth-data"])
        self.assertIn("1 dropped beyond it", out)

    def test_counts_sum_to_births_inside(self):
        out = run([self.path, "--col", "x=x", "--col", "z=z", "--col", "t=t",
                   "--i-have-birth-data"])
        counts = eval(out.split("per-bin counts:")[1].split("\n")[0])
        self.assertEqual(sum(counts), 451)
        self.assertEqual(counts, [141, 56, 46, 101, 87, 20])


class MovingAxis(unittest.TestCase):
    def test_axis_drift_is_applied(self):
        # A pore sitting exactly on the axis at t=1ms is at x = 250um + 36um.
        p = write_csv([(286e-6, 250e-6, 1e-3)], ["x", "z", "t"])
        try:
            out = run([p, "--col", "x=x", "--col", "z=z", "--col", "t=t",
                       "--i-have-birth-data"])
            self.assertIn("[1, 0, 0, 0, 0, 0]", out)
        finally:
            os.unlink(p)

    def test_ignoring_drift_puts_it_in_the_wrong_bin(self):
        # Same pore, axis held still: it lands 36 um out, in the 20-40 bin.
        p = write_csv([(286e-6, 250e-6, 1e-3)], ["x", "z", "t"])
        try:
            out = run([p, "--col", "x=x", "--col", "z=z", "--col", "t=t",
                       "--axis-velocity", "0", "--i-have-birth-data"])
            self.assertIn("[0, 1, 0, 0, 0, 0]", out)
        finally:
            os.unlink(p)

    def test_moving_axis_without_time_is_refused(self):
        p = write_csv([(286e-6, 250e-6)], ["x", "z"])
        try:
            with self.assertRaises(SystemExit):
                run([p, "--col", "x=x", "--col", "z=z", "--i-have-birth-data"])
        finally:
            os.unlink(p)


class RadiusColumn(unittest.TestCase):
    def test_radius_given_directly(self):
        p = write_csv([(1e-5,), (3e-5,), (9e-5,)], ["r"])
        try:
            out = run([p, "--col", "r=r", "--i-have-birth-data"])
            self.assertIn("[1, 1, 0, 0, 1, 0]", out)
        finally:
            os.unlink(p)

    def test_micron_scaling(self):
        p = write_csv([(10,), (30,), (90,)], ["r"])
        try:
            out = run([p, "--col", "r=r", "--scale-um", "--i-have-birth-data"])
            self.assertIn("[1, 1, 0, 0, 1, 0]", out)
        finally:
            os.unlink(p)


class Provenance(unittest.TestCase):
    def test_must_declare_which_kind_of_data(self):
        p = write_csv([(1e-5,)], ["r"])
        try:
            with self.assertRaises(SystemExit):
                M.main([p, "--col", "r=r"])
        finally:
            os.unlink(p)

    def test_final_positions_emits_the_warning(self):
        p = write_csv([(1e-5,)], ["r"])
        try:
            out = run([p, "--col", "r=r", "--final-positions"])
            self.assertIn("WARNING", out)
            self.assertIn("circular", out)
        finally:
            os.unlink(p)

    def test_birth_data_emits_no_warning(self):
        p = write_csv([(1e-5,)], ["r"])
        try:
            out = run([p, "--col", "r=r", "--i-have-birth-data"])
            self.assertNotIn("WARNING", out)
        finally:
            os.unlink(p)


class SizeDistribution(unittest.TestCase):
    def test_grid_is_uniform(self):
        """The trapezoid bug: `general` reads the table as a density."""
        p = write_csv([(1e-5, 1.1e-5), (1e-5, 1.5e-5), (1e-5, 3.1e-5)],
                      ["r", "d"])
        try:
            out = run([p, "--col", "r=r", "--col", "d=d", "--i-have-birth-data"])
            xs = [float(l.strip().split()[0].lstrip("("))
                  for l in out.splitlines()
                  if l.strip().startswith("(") and "e-" in l and "0." in l]
            steps = np.diff(sorted(xs))
            np.testing.assert_allclose(steps, M.DIAM_STEP, atol=1e-12)
        finally:
            os.unlink(p)

    def test_absent_diameter_column_says_so(self):
        p = write_csv([(1e-5,)], ["r"])
        try:
            out = run([p, "--col", "r=r", "--i-have-birth-data"])
            self.assertIn("sizeDistribution left unchanged", out)
        finally:
            os.unlink(p)


class Rate(unittest.TestCase):
    def test_rate_is_per_liquid_volume_per_second(self):
        rows = [(1e-5,)] * 100
        p = write_csv(rows, ["r"])
        try:
            out = run([p, "--col", "r=r", "--melt-volume", "1e-11",
                       "--window", "1e-3", "--i-have-birth-data"])
            line = next(l for l in out.splitlines() if l.startswith("rate"))
            self.assertAlmostEqual(float(line.split()[1].rstrip(";")), 1e16, delta=1e13)
        finally:
            os.unlink(p)

    def test_no_rate_without_both_inputs(self):
        p = write_csv([(1e-5,)], ["r"])
        try:
            out = run([p, "--col", "r=r", "--melt-volume", "1e-11",
                       "--i-have-birth-data"])
            self.assertFalse(any(l.startswith("rate") for l in out.splitlines()))
        finally:
            os.unlink(p)


class BadInput(unittest.TestCase):
    def test_missing_column_names_the_alternatives(self):
        p = write_csv([(1e-5,)], ["radius"])
        try:
            with self.assertRaises(SystemExit) as cm:
                M.main([p, "--col", "r=nope", "--i-have-birth-data"])
            self.assertIn("radius", str(cm.exception))
        finally:
            os.unlink(p)

    def test_all_births_outside_outer_edge_is_an_error(self):
        p = write_csv([(5e-4,)], ["r"])
        try:
            with self.assertRaises(SystemExit):
                M.main([p, "--col", "r=r", "--i-have-birth-data"])
        finally:
            os.unlink(p)

    def test_runs_as_a_subprocess(self):
        p = write_csv([(1e-5,), (9e-5,)], ["r"])
        try:
            r = subprocess.run([sys.executable, SCRIPT, p, "--col", "r=r",
                                "--i-have-birth-data"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("radialSeeding", r.stdout)
        finally:
            os.unlink(p)


if __name__ == "__main__":
    unittest.main(verbosity=2)
