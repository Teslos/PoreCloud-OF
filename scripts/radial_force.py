"""Is the inward pore motion driven by JxB, or is it an artefact of the moving laser?

Two separate checks, because the first analysis conflated them.

(1) THE METRIC. Radial displacement was measured against the INSTANTANEOUS laser
    position, and the laser advances 36 um during the 0.975 ms window. A pore that
    never moves still shows r shrinking by up to 36 um if it sits ahead of the
    beam -- comparable to the -34 um that was reported as inward transport, and
    identical across cases, which would explain why every configuration looked
    alike. So recompute against a FIXED origin (the laser position at the pore's
    birth) and against the lab frame, and report the difference.

(2) THE FORCE. Displacement mixes advection, buoyancy and JxB. The force itself
    does not. Decompose the stored exclusion force into horizontal radial and
    azimuthal components about the keyhole axis:

        F_r   = F_h . r_hat      (< 0 = pushing the pore toward the cavity)
        F_th  = F_h . th_hat     (swirl)

    For B along y, JxB = B_y(-J_z, 0, J_x): the horizontal current rotated 90
    degrees about the vertical. A purely RADIAL current therefore gives a purely
    AZIMUTHAL force and no radial transport at all; radial force requires an
    azimuthal current, which an axisymmetric keyhole does not have. Whether the
    scanning asymmetry supplies one is an empirical question -- hence this script.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.expanduser("~"))
import analyse_cloud as ac


def decompose(d):
    """Radial/azimuthal split of the exclusion force about the keyhole axis."""
    rx = d["Cx_m"] - ac.laser_x(d["Time"])
    rz = d["Cz_m"] - ac.ZAXIS
    r = np.hypot(rx, rz)
    ok = r > 1e-9
    rhx, rhz = rx[ok] / r[ok], rz[ok] / r[ok]
    Fx, Fz = d["Fx_N"][ok], d["Fz_N"][ok]
    Fr = Fx * rhx + Fz * rhz                 # + outward, - inward
    Fth = -Fx * rhz + Fz * rhx               # azimuthal
    return Fr, Fth, r[ok]


def main(cases):
    print("=" * 108)
    print("(1) HOW MUCH OF THE INWARD MOTION IS THE LASER MOVING?")
    print("=" * 108)
    print("%-26s %10s %10s %10s %10s %10s" % (
        "case", "dr_moving", "dr_fixed", "laser part", "dx_lab", "dz_lab"))
    print("%-26s %10s %10s %10s %10s %10s" % ("", "um", "um", "um", "um", "um"))
    store = {}
    for c in cases:
        d = ac.load(os.path.expanduser("~/results/%s/poreCloud_pores.csv" % c))
        if d is None:
            continue
        p = ac.per_pore(d)
        # against the MOVING laser (what was reported before)
        r0m = np.hypot(p["Cx_m_0"] - ac.laser_x(p["t0"]), p["Cz_m_0"] - ac.ZAXIS)
        r1m = np.hypot(p["Cx_m_1"] - ac.laser_x(p["t1"]), p["Cz_m_1"] - ac.ZAXIS)
        # against a FIXED origin: the laser position when that pore was born
        r0f = r0m
        r1f = np.hypot(p["Cx_m_1"] - ac.laser_x(p["t0"]), p["Cz_m_1"] - ac.ZAXIS)
        drm, drf = (r1m - r0m) * 1e6, (r1f - r0f) * 1e6
        dx = (p["Cx_m_1"] - p["Cx_m_0"]) * 1e6
        dz = (p["Cz_m_1"] - p["Cz_m_0"]) * 1e6
        print("%-26s %+10.2f %+10.2f %+10.2f %+10.2f %+10.2f" % (
            c.replace("poreCloud-", ""), np.median(drm), np.median(drf),
            np.median(drm) - np.median(drf), np.median(dx), np.median(dz)))
        store[c] = d

    print()
    print("=" * 108)
    print("(2) THE JxB FORCE ITSELF: radial vs azimuthal about the keyhole axis")
    print("=" * 108)
    print("%-26s %12s %8s %12s %8s %9s" % (
        "case", "med F_r", "inward", "med |F_th|", "azim>rad", "|F_th|/|F_r|"))
    print("%-26s %12s %8s %12s %8s %9s" % ("", "N", "%", "N", "%", ""))
    for c, d in store.items():
        Fr, Fth, r = decompose(d)
        med_fr = np.median(Fr)
        inward = 100.0 * (Fr < 0).mean()
        azim = 100.0 * (np.abs(Fth) > np.abs(Fr)).mean()
        ratio = np.median(np.abs(Fth)) / np.median(np.abs(Fr))
        print("%-26s %+12.3e %8.1f %12.3e %8.1f %9.2f" % (
            c.replace("poreCloud-", ""), med_fr, inward,
            np.median(np.abs(Fth)), azim, ratio))

    print()
    print("  F_r < 0 means the force points TOWARD the keyhole axis.")
    print("  50% inward = no radial preference; the force is purely swirl on average.")
    print("  |F_th|/|F_r| >> 1 means the force is predominantly azimuthal (stirring),")
    print("  which is what a vertical B gives when the thermoelectric current is radial.")


if __name__ == "__main__":
    main(["poreCloud-ac-by-400Hz-fix", "poreCloud-rmf-xz-0.2T", "poreCloud-rmf-yz-0.2T",
          "poreCloud-dc-bx-0.2T", "poreCloud-azimuthal-1T",
          "poreCloud-rmf-yz-0.2T-seeded"])
