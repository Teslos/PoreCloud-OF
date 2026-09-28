"""Score a PoreCloud Lagrangian run, and in particular whether pores move TOWARD
the keyhole.

Motivated by the experimental result for exp `bz` (= sim `by`, vertical): pores are
driven laterally INTO the keyhole cavity rather than vented upward
(esrf-pore-fate-tracking). For a vertical field that is the only thing the
exclusion force can do, because (JxB)_y = J_z B_x - J_x B_z vanishes identically
when B is along y -- so the field cannot lift or sink a pore at all, and any
vertical motion is advection plus buoyancy.

The keyhole axis is vertical, under the (moving) laser spot, so the experimental
claim becomes a statement about the HORIZONTAL radial coordinate

    r(t) = |(Cx - x_laser(t), Cz - z_laser)|

and the test is the sign of dr from birth to last sighting. dr < 0 means driven
inward, toward the cavity -- the simulated counterpart of "absorbed at the keyhole".

Reports, per case: net displacement by component, the radial test, pore speed in
the experiment's units (um/ms), and the vertical force budget, which for a
vertical field contains no electromagnetic term at all.
"""
import csv
import os
import sys

import numpy as np

# laser programme, shared by every Lagrangian case (constant/timeVsLaserPosition):
# (0, (250e-6 700e-6 250e-6)) -> (2.5e-3, (340e-6 700e-6 250e-6)) = 36 mm/s along x
X0, XEND, TEND = 250e-6, 340e-6, 2.5e-3
ZAXIS = 250e-6
YSURF = 700e-6


def laser_x(t):
    return X0 + (XEND - X0) * np.clip(t / TEND, 0, 1)


def load(path):
    rows = list(csv.DictReader(open(path)))
    if not rows:
        return None
    d = {}
    for k in rows[0]:
        try:
            d[k] = np.array([float(r[k]) for r in rows])
        except ValueError:
            d[k] = np.array([r[k] for r in rows])
    return d


def per_pore(d):
    """Collapse the time series to one record per pore: birth and last sighting."""
    pid = d["PoreID"].astype(int)
    t = d["Time"]
    order = np.lexsort((t, pid))
    pid, t = pid[order], t[order]
    idx = {k: d[k][order] for k in ("Cx_m", "Cy_m", "Cz_m", "d_m") if k in d}
    first = np.concatenate(([True], pid[1:] != pid[:-1]))
    last = np.concatenate((pid[1:] != pid[:-1], [True]))
    out = {"t0": t[first], "t1": t[last], "n": int(first.sum())}
    for k in idx:
        out[k + "_0"] = idx[k][first]
        out[k + "_1"] = idx[k][last]
    return out


def q(v, name, unit="um", scale=1e6):
    v = np.asarray(v) * scale
    return "%-22s median %+8.2f   IQR %+8.2f .. %+8.2f   mean %+8.2f  [%s]" % (
        name, np.median(v), *np.percentile(v, [25, 75]), v.mean(), unit)


def report(path, label):
    d = load(path)
    if d is None:
        print("%-28s no rows" % label)
        return
    p = per_pore(d)
    print("=" * 94)
    print("%s" % label)
    print("=" * 94)
    span = (d["Time"].max() - d["Time"].min()) * 1e3
    print("  %d pores, %d samples, t = %.3f .. %.3f ms (span %.3f ms)"
          % (p["n"], len(d["Time"]), d["Time"].min() * 1e3, d["Time"].max() * 1e3, span))

    # --- force components: the F.B = 0 test ------------------------------------
    if "Fy_N" in d:
        print("\n  exclusion force, max|F| per component  (B along y => F_y must vanish)")
        print("    |Fx| %.3e   |Fy| %.3e   |Fz| %.3e  N"
              % (np.abs(d["Fx_N"]).max(), np.abs(d["Fy_N"]).max(),
                 np.abs(d["Fz_N"]).max()))
        if "Fbuoyy" in d:
            fb = np.abs(d["Fbuoyy"])
            fe = np.sqrt(d["Fx_N"] ** 2 + d["Fy_N"] ** 2 + d["Fz_N"] ** 2)
            nz = fb > 0
            print("    |F_excl|/|F_buoy| median %.2f (n=%d)"
                  % (np.median(fe[nz] / fb[nz]), nz.sum()))
            print("    VERTICAL budget: F_y(EM) = %.1e, so the only vertical force is"
                  " buoyancy, median |F_buoy_y| = %.3e N"
                  % (np.abs(d["Fy_N"]).max(), np.median(fb[nz])))

    # --- net displacement -----------------------------------------------------
    dy = p["Cy_m_1"] - p["Cy_m_0"]
    dx = p["Cx_m_1"] - p["Cx_m_0"]
    dz = p["Cz_m_1"] - p["Cz_m_0"]
    print("\n  net displacement, birth -> last sighting")
    print("    " + q(dy, "dy  (vertical, +up)"))
    print("    " + q(dx, "dx  (scan)"))
    print("    " + q(dz, "dz  (transverse)"))

    # --- the radial test: driven toward the keyhole? ---------------------------
    r0 = np.hypot(p["Cx_m_0"] - laser_x(p["t0"]), p["Cz_m_0"] - ZAXIS)
    r1 = np.hypot(p["Cx_m_1"] - laser_x(p["t1"]), p["Cz_m_1"] - ZAXIS)
    dr = r1 - r0
    print("\n  RADIAL test vs the keyhole axis  (dr < 0 = driven INWARD, toward the cavity)")
    print("    " + q(dr, "dr"))
    print("    inward fraction %.1f%%  (%d of %d pores)"
          % (100.0 * (dr < 0).mean(), int((dr < 0).sum()), p["n"]))
    print("    median r: birth %.1f um -> last %.1f um" % (np.median(r0) * 1e6,
                                                          np.median(r1) * 1e6))

    # --- speed, in the experiment's units -------------------------------------
    if "Ux" in d:
        sp = np.sqrt(d["Ux"] ** 2 + d["Uy"] ** 2 + d["Uz"] ** 2) * 1e3  # m/s -> um/ms
        print("\n  pore speed (experiment reports um/ms; exp bz 77.3, no-field 51.6)")
        print("    median %.1f  IQR %.1f .. %.1f  um/ms"
              % (np.median(sp), *np.percentile(sp, [25, 75])))
        net = np.hypot(np.hypot(dx, dy), dz) / ((p["t1"] - p["t0"]) + 1e-12) * 1e3
        ok = (p["t1"] - p["t0"]) > 1e-5
        if ok.sum():
            print("    net transport speed |d|/dt: median %.1f um/ms (n=%d)"
                  % (np.median(net[ok]), ok.sum()))
    if "Uslipy" in d:
        print("    slip |Uslip| median %.2f um/ms"
              % (np.median(np.sqrt(d["Uslipx"] ** 2 + d["Uslipy"] ** 2
                                   + d["Uslipz"] ** 2)) * 1e3))

    # --- sizes, for the published-population comparison -----------------------
    if "d_m_0" in p:
        dd = p["d_m_0"] * 1e6
        print("\n  pore diameter: median %.1f um, range %.1f .. %.1f"
              % (np.median(dd), dd.min(), dd.max()))
    print()


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        path = os.path.join(os.path.expanduser(arg), "poreCloud_pores.csv")
        if os.path.exists(path):
            report(path, os.path.basename(arg.rstrip("/")))
        else:
            print("%-28s no poreCloud_pores.csv" % os.path.basename(arg.rstrip("/")))
