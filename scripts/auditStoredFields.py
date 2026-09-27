"""Decide from the STORED FORCE FIELDS what each case actually ran with.

Needed because in three of four case families every directory carries an
identical MHD dict -- same B0, same frequency, same rotationPlane, and
`active true` even in no-mhd -- so the dicts cannot say which field a run used.
The force fields are data written by the run itself, so they can.

For each case, over the liquid region: mean |F| and the share of force-squared
along each axis. A rotating field confined to a plane drives a force with a
characteristic component split, so two cases that really ran different planes
must differ here. Cases that agree to rounding ran the same configuration.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flowStructureSurvey import read_field, set_mesh

fam = os.path.expanduser(sys.argv[1])
t = sys.argv[2] if len(sys.argv) > 2 else "0.0005"
ctl = os.path.join(fam, "no-mhd")

cc = [read_field(os.path.join(ctl, "0", "C" + c), 1) for c in "zyx"]
coords = [np.unique(v) for v in cc]
set_mesh(len(coords[2]), len(coords[1]), len(coords[0]))

print("family %s   t=%s" % (fam, t))
print("%-16s %12s %8s %8s %8s %10s" % (
    "case", "mean|F|", "fx2", "fy2", "fz2", "sha1-ish"))
ref = {}
for c in sorted(d for d in os.listdir(fam) if os.path.isdir(os.path.join(fam, d))):
    p = os.path.join(fam, c, t, "LorentzForce")
    F = read_field(p, 3)
    al = read_field(os.path.join(fam, c, t, "alpha.metal"), 1)
    T = read_field(os.path.join(fam, c, t, "T"), 1)
    TL = read_field(os.path.join(fam, c, t, "TLiquidus"), 1)
    if F is None or al is None:
        print("%-16s %12s" % (c, "no field file"))
        continue
    liq = (al > 0.5) & (T > (TL if TL is not None else 867.0))
    if liq.sum() < 100:
        print("%-16s %12s" % (c, "no liquid"))
        continue
    Fl = F[liq]
    mag = np.linalg.norm(Fl, axis=1)
    s2 = (Fl ** 2).sum(0)
    tot = s2.sum()
    # cheap content fingerprint: two runs with the same field give the same digits
    fp = "%.6e" % np.abs(Fl).sum()
    print("%-16s %12.4g %8.3f %8.3f %8.3f %10s" % (
        c, mag.mean(), *(s2 / tot if tot > 0 else [np.nan] * 3), fp[-9:]))
    ref[c] = fp

dups = {}
for c, fp in ref.items():
    dups.setdefault(fp, []).append(c)
same = [v for v in dups.values() if len(v) > 1]
if same:
    print("\nIDENTICAL force fields -- these ran the SAME configuration:")
    for g in same:
        print("   " + " == ".join(sorted(g)))
else:
    print("\nAll cases have distinct force fields: each ran its own configuration.")
