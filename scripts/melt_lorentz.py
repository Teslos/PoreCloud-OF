"""Does the Lorentz force on the MELT drive flow toward the keyhole?

The previous check decomposed the exclusion force on the PORES and found it
predominantly azimuthal, with a median radial component pointing outward -- yet
pores move inward. So the inward transport is not the pore force; it must be
advection. That leaves the other JxB term: the Lorentz force in the momentum
equation, acting on the liquid, which reorganises the flow that then carries pores.

These are different quantities and the project has not separated them:

    exclusion force on a bubble   F = -(3/4) V (JxB)     -> Lagrangian, per pore
    body force on the melt            JxB                -> Eulerian, per volume

This script takes the stored Eulerian `LorentzForce` and `U` fields and reports,
over the liquid region and about the vertical keyhole axis under the laser:

    F_r    radial component of the melt body force   (< 0 = pushing melt inward)
    U_r    radial component of the melt velocity     (< 0 = flow running inward)
    F_th   azimuthal component of the body force     (swirl drive)

If a vertical field drives inward MELT flow, F_r should be negative for the AC-by
case and that should show up in U_r relative to the other configurations.
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser("~"))
from flowStructureSurvey import grid, read_field, set_mesh

GEOM_T = "0.001"
ZAXIS, YSURF = 250e-6, 700e-6
X0, XEND, TEND = 250e-6, 340e-6, 2.5e-3


def laser_x(t):
    return X0 + (XEND - X0) * min(max(t / TEND, 0.0), 1.0)


def case_times(case):
    return sorted((d for d in os.listdir(case) if re.fullmatch(r"[0-9.eE+-]+", d)),
                  key=float)


def setup(case):
    bmd = os.path.join(case, "system", "blockMeshDict")
    m = re.search(r"hex\s*\([^)]*\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", open(bmd).read())
    nx, ny, nz = (int(g) for g in m.groups())
    set_mesh(nx, ny, nz)
    cc = [read_field(os.path.join(case, GEOM_T, "C" + c), 1) for c in "xyz"]
    if any(v is None for v in cc):
        return None
    g = [grid(v) for v in cc]
    return g[0][0, 0, :], g[1][0, :, 0], g[2][:, 0, 0]


def score(case, t):
    """Radial/azimuthal split of the melt body force and velocity, in the liquid."""
    F = read_field(os.path.join(case, t, "LorentzForce"), 3)
    U = read_field(os.path.join(case, t, "U"), 3)
    al = read_field(os.path.join(case, t, "alpha.metal"), 1)
    T = read_field(os.path.join(case, t, "T"), 1)
    TL = read_field(os.path.join(case, t, "TLiquidus"), 1)
    if F is None or U is None or al is None:
        return None
    liq = (al > 0.5) & (T > (TL if TL is not None else 867.0))
    if liq.sum() < 500:
        return None

    xc, yc, zc = setup(case)
    X = np.broadcast_to(xc.reshape(1, 1, -1), (len(zc), len(yc), len(xc))).ravel()
    Z = np.broadcast_to(zc.reshape(-1, 1, 1), (len(zc), len(yc), len(xc))).ravel()
    rx, rz = X[liq] - laser_x(float(t)), Z[liq] - ZAXIS
    r = np.hypot(rx, rz)
    ok = r > 1e-9
    rhx, rhz = rx[ok] / r[ok], rz[ok] / r[ok]

    Fl, Ul = F[liq][ok], U[liq][ok]
    Fr = Fl[:, 0] * rhx + Fl[:, 2] * rhz
    Fth = -Fl[:, 0] * rhz + Fl[:, 2] * rhx
    Ur = Ul[:, 0] * rhx + Ul[:, 2] * rhz
    return dict(n=int(ok.sum()), Fr=np.median(Fr), Fr_in=100.0 * (Fr < 0).mean(),
                Fth=np.median(np.abs(Fth)), Fabs=np.median(np.abs(Fr)),
                Ur=np.median(Ur), Ur_in=100.0 * (Ur < 0).mean())


def main(cases):
    print("%-26s %5s %12s %8s %12s %8s %10s %8s" % (
        "case", "t/ms", "med F_r", "F inward", "med|F_th|", "|Fth/Fr|", "med U_r", "U inward"))
    print("%-26s %5s %12s %8s %12s %8s %10s %8s" % (
        "", "", "N/m3", "%", "N/m3", "", "m/s", "%"))
    for c in cases:
        case = os.path.expanduser("~/results/%s" % c)
        ts = [t for t in case_times(case) if float(t) > 0]
        if not ts:
            print("%-26s no reconstructed times" % c.replace("poreCloud-", ""))
            continue
        for t in [GEOM_T]:
            s = score(case, t)
            if s is None:
                print("%-26s %5s  unusable" % (c.replace("poreCloud-", ""), t))
                continue
            print("%-26s %5.2f %+12.3e %8.1f %12.3e %8.2f %+10.4f %8.1f" % (
                c.replace("poreCloud-", ""), float(t) * 1e3, s["Fr"], s["Fr_in"],
                s["Fth"], s["Fth"] / max(s["Fabs"], 1e-30), s["Ur"], s["Ur_in"]))
    print()
    print("  F_r < 0 = body force pushes melt TOWARD the keyhole axis.")
    print("  U_r < 0 = melt is actually flowing inward there.")
    print("  50% = no radial preference. |Fth/Fr| >> 1 = the drive is swirl, not radial.")


if __name__ == "__main__":
    main(["poreCloud-ac-by-400Hz-fix", "poreCloud-rmf-xz-0.2T", "poreCloud-rmf-yz-0.2T",
          "poreCloud-dc-bx-0.2T", "poreCloud-azimuthal-1T"])
