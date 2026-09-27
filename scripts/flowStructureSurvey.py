"""Score CHANGE IN FLOW STRUCTURE between a field-on case and its field-off control.

Every internal sweep so far has scored the FORCE DIRECTION on pores --
cos(<F>,+y), or the fraction of liquid in which JxB points up. Two separate
results have now undermined that instrument:

  * melt-pool-vortex-circulation: pores escape by riding a vortex back to the
    surface, not by being pushed up, so alignment of force with +y mis-scores
    lateral transport;
  * the F.B = 0 identity: a field exactly along an axis exerts no force along
    that axis, so a vertical field is UNSCOREABLE by construction -- the survey
    records dc-by-0.2T as exactly 0.0, which is an artefact of the metric, not
    a measurement.

Fan's demonstrated mechanism is reorganisation of the melt-pool FLOW. This
script scores that instead, from the stored velocity fields of runs already
performed. Nothing is re-simulated.

Metrics, all on the liquid region only (alpha.metal > 0.5 and T > TLiquidus):

  speed      <|u|>                     flow intensity: does the field drive or brake
  omega      <|curl u|>                rotational intensity
  fx,fy,fz   <w_i^2>/<|w|^2>           VORTEX AXIS FINGERPRINT. w_z is circulation
                                       in the scan-vertical plane, i.e. the plane a
                                       pore must circulate in to reach the surface.
  up_frac    V(u_y>0)/V                fraction of liquid FLOWING up. Distinct from
  up_flux    <max(u_y,0)>              the old metric, which used force direction.
  Qpos       V(Q>0)/V                  volume occupied by vortex cores
  D_cos      1 - <cos(u_on,u_off)>     direct structural distance from the control,
  D_L2       |u_on-u_off|/|u_off|      on cells liquid in BOTH cases

Vorticity and Q are evaluated on the liquid mask ERODED by one cell: the
solid/liquid boundary carries an artificial shear layer (Darcy damping drives u
to zero inside the solid) and un-eroded gradients there swamp the interior.

D_cos and D_L2 cannot by themselves distinguish a field effect from chaotic
divergence of two nearby trajectories. Pass --chaos to bound that: it scores
0.1 T against 0.2 T of the SAME rotation plane, which is a field change of the
same kind, so a D of comparable size means the metric is saturated.

    python3 flow_structure.py --family ~/results_matched_0.2T_f4000hz
"""

import argparse
import os
import re
import sys

import numpy as np

# Set per family by set_mesh() from the cell-centre fields; never hardcoded,
# because reshaping to the wrong dims yields plausible nonsense, not an error.
NX = NY = NZ = NCELLS = None


def set_mesh(nx, ny, nz):
    global NX, NY, NZ, NCELLS
    NX, NY, NZ = nx, ny, nz
    NCELLS = nx * ny * nz
    return NCELLS


# ---------------------------------------------------------------- IO
def read_field(path, ncomp, ncells=None):
    """OpenFOAM volScalar/volVectorField -> (ncells,) or (ncells,3).

    Handles both `writeFormat ascii` and `writeFormat binary`: the header is text
    in either case, but a binary payload is raw little-endian float64 between the
    parentheses. Everything is read as bytes so a binary body cannot crash the
    text decode. Uniform fields are broadcast.
    """
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        raw = f.read()
    i = raw.find(b"internalField")
    if i < 0:
        return None
    n = ncells if ncells is not None else NCELLS

    m = re.match(rb"internalField\s+uniform\s+(.+?);", raw[i:i + 200], re.S)
    if m:
        v = np.fromstring(m.group(1).replace(b"(", b" ").replace(b")", b" ")
                          .decode("ascii"), sep=" ")
        if n is None:
            return None
        return np.tile(v, (n, 1)) if ncomp == 3 else np.full(n, float(v[0]))

    lst = raw.index(b"List<", i)
    cnt = int(re.search(rb"(\d+)", raw[raw.index(b">", lst):lst + 120]).group(1))
    if n is None:
        n = cnt
    elif cnt != n:
        raise ValueError("%s: declares %d cells, mesh has %d" % (path, cnt, n))
    j = raw.index(b"(", lst)

    header = raw[:i]
    binary = b"binary" in (header.split(b"format", 1)[1][:40]
                           if b"format" in header else b"")
    if binary:
        a = np.frombuffer(raw, dtype="<f8", count=n * ncomp, offset=j + 1)
    else:
        k = raw.index(b"\n)", j)
        body = raw[j + 1:k]
        if ncomp == 3:
            body = body.replace(b"(", b" ").replace(b")", b" ")
        a = np.fromstring(body.decode("ascii", "replace"), sep=" ")
    if a.size != n * ncomp:
        raise ValueError("%s: got %d values, want %d" % (path, a.size, n * ncomp))
    return a.reshape(n, 3) if ncomp == 3 else np.asarray(a, dtype=float)


def grid(a):
    """flat OpenFOAM cell order -> (nz, ny, nx[, 3]); i fastest, then j, then k."""
    return a.reshape(NZ, NY, NX) if a.ndim == 1 else a.reshape(NZ, NY, NX, 3)


def times(case, skip_frac=0.2):
    """Sorted write times, dropping the startup transient."""
    ts = sorted((d for d in os.listdir(case) if re.fullmatch(r"[0-9.eE+-]+", d)),
                key=float)
    ts = [t for t in ts if float(t) > 0]
    return ts[int(len(ts) * skip_frac):]


# ---------------------------------------------------------------- geometry
def erode(mask):
    """Drop cells with a non-liquid 6-neighbour: the interface shear layer is
    numerical (Darcy damping), and its gradients dwarf the interior flow."""
    e = mask.copy()
    for ax in (0, 1, 2):
        for sh in (1, -1):
            e &= np.roll(mask, sh, axis=ax)
    e[0] = e[-1] = False
    e[:, 0] = e[:, -1] = False
    e[:, :, 0] = e[:, :, -1] = False
    return e


def grad_tensor(u, zc, yc, xc):
    """J[i][j] = du_i/dx_j with x=(x,y,z); handles the graded spacing."""
    d = [np.gradient(u[..., i], zc, yc, xc) for i in range(3)]  # [dz,dy,dx]
    return [[d[i][2], d[i][1], d[i][0]] for i in range(3)]      # -> [dx,dy,dz]


# ---------------------------------------------------------------- metrics
def snapshot(case, t, coords, vol):
    zc, yc, xc = coords
    U = read_field(os.path.join(case, t, "U"), 3)
    al = read_field(os.path.join(case, t, "alpha.metal"), 1)
    T = read_field(os.path.join(case, t, "T"), 1)
    TL = read_field(os.path.join(case, t, "TLiquidus"), 1)
    if U is None or al is None or T is None:
        return None
    if TL is None:
        TL = np.full(NCELLS, 867.0)          # AlSi10Mg liquidus, fallback only

    liq = grid((al > 0.5) & (T > TL))
    if liq.sum() < 500:
        return None
    u, V = grid(U), grid(vol)
    w = V[liq].sum()

    core = erode(liq)
    out = {"t": float(t), "n_liq": int(liq.sum()), "n_core": int(core.sum()),
           "V_liq_mm3": w * 1e9}

    # --- intensity and vertical transport (whole liquid) ---
    sp = np.linalg.norm(u, axis=-1)
    out["speed"] = float((sp[liq] * V[liq]).sum() / w)
    uy = u[..., 1]
    out["up_frac"] = float(V[liq & (uy > 0)].sum() / w)
    out["up_flux"] = float((np.maximum(uy[liq], 0) * V[liq]).sum() / w)

    if out["n_core"] < 200:
        return out

    # --- rotational structure (eroded core only) ---
    J = grad_tensor(u, zc, yc, xc)
    wx = J[2][1] - J[1][2]
    wy = J[0][2] - J[2][0]
    wz = J[1][0] - J[0][1]
    Vc, wc = V[core], V[core].sum()
    w2 = wx[core] ** 2 + wy[core] ** 2 + wz[core] ** 2
    out["omega"] = float((np.sqrt(w2) * Vc).sum() / wc)
    tot = (w2 * Vc).sum()
    for nm, comp in (("fx", wx), ("fy", wy), ("fz", wz)):
        out[nm] = float(((comp[core] ** 2) * Vc).sum() / tot) if tot > 0 else np.nan

    S2 = O2 = 0.0
    for i in range(3):
        for j in range(3):
            s = 0.5 * (J[i][j] + J[j][i])
            o = 0.5 * (J[i][j] - J[j][i])
            S2 = S2 + s[core] ** 2
            O2 = O2 + o[core] ** 2
    Q = 0.5 * (O2 - S2)
    out["Qpos"] = float(Vc[Q > 0].sum() / wc)
    return out


def compare(case_a, case_b, t, coords, vol):
    """Structural distance between two cases at the same time, same mesh."""
    fa = [read_field(os.path.join(case_a, t, f), n)
          for f, n in (("U", 3), ("alpha.metal", 1), ("T", 1), ("TLiquidus", 1))]
    fb = [read_field(os.path.join(case_b, t, f), n)
          for f, n in (("U", 3), ("alpha.metal", 1), ("T", 1), ("TLiquidus", 1))]
    if fa[0] is None or fb[0] is None:
        return None

    def mask(f):
        TL = f[3] if f[3] is not None else np.full(NCELLS, 867.0)
        return (f[1] > 0.5) & (f[2] > TL)

    m = mask(fa) & mask(fb)
    if m.sum() < 500:
        return None
    ua, ub, V = fa[0][m], fb[0][m], vol[m]
    na, nb = np.linalg.norm(ua, axis=1), np.linalg.norm(ub, axis=1)
    ok = (na > 1e-9) & (nb > 1e-9)
    cos = ((ua[ok] * ub[ok]).sum(1) / (na[ok] * nb[ok]))
    d_l2 = (np.sqrt((V[:, None] * (ua - ub) ** 2).sum())
            / np.sqrt((V[:, None] * ub ** 2).sum()))
    return {"t": float(t), "n_both": int(m.sum()),
            "D_cos": float(1.0 - (cos * V[ok]).sum() / V[ok].sum()),
            "D_L2": float(d_l2)}


# ---------------------------------------------------------------- driver
def med_iqr(rows, key):
    v = np.array([r[key] for r in rows if key in r and np.isfinite(r[key])])
    if not v.size:
        return np.nan, np.nan, np.nan
    return np.median(v), np.percentile(v, 25), np.percentile(v, 75)


def paired(rows_a, rows_b, key):
    """Wilcoxon signed-rank on time-matched snapshot pairs. Returns (ratio, p, n)."""
    from scipy.stats import wilcoxon
    ta = {r["t"]: r.get(key) for r in rows_a}
    tb = {r["t"]: r.get(key) for r in rows_b}
    pairs = [(ta[t], tb[t]) for t in sorted(set(ta) & set(tb))
             if ta.get(t) is not None and tb.get(t) is not None
             and np.isfinite(ta[t]) and np.isfinite(tb[t])]
    if len(pairs) < 6:
        return np.nan, np.nan, len(pairs)
    x = np.array([q[0] for q in pairs]); y = np.array([q[1] for q in pairs])
    ratio = np.median(x) / np.median(y) if np.median(y) else np.nan
    if np.allclose(x, y):
        return ratio, 1.0, len(pairs)
    return ratio, float(wilcoxon(x, y).pvalue), len(pairs)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--family", required=True)
    p.add_argument("--control", default="no-mhd")
    p.add_argument("--max-snaps", type=int, default=20)
    p.add_argument("--chaos", action="store_true",
                   help="also score 0.1 T vs 0.2 T of the same plane, to bound "
                        "how much of D is chaotic divergence")
    a = p.parse_args()

    fam = os.path.expanduser(a.family)
    ctl = os.path.join(fam, a.control)
    cases = sorted(d for d in os.listdir(fam)
                   if os.path.isdir(os.path.join(fam, d))
                   and (d.startswith(("rmf", "dc")) or d == a.control))

    bmd = os.path.join(ctl, "system", "blockMeshDict")
    if not os.path.exists(bmd):
        bmd = os.path.join(ctl, "constant", "polyMesh", "blockMeshDict")
    m = re.search(r"hex\s*\([^)]*\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)",
                  open(bmd).read())
    if not m:
        sys.exit("cannot read (nx ny nz) from %s" % bmd)
    nx, ny, nz = (int(g) for g in m.groups())
    n = set_mesh(nx, ny, nz)

    cc = [read_field(os.path.join(ctl, "0", "C" + c), 1) for c in "xyz"]
    if any(v is None for v in cc):
        sys.exit("no cell-centre fields -- run prep_geom.sh for this family")
    if cc[0].size != n:
        sys.exit("blockMeshDict says %d cells, Cx has %d" % (n, cc[0].size))

    # one line along each axis; grading round-off within a column is not a
    # different coordinate, so check it is negligible rather than exact
    gx, gy, gz = (grid(v) for v in cc)
    xc, yc, zc = gx[0, 0, :], gy[0, :, 0], gz[:, 0, 0]
    for nm, g, ref, ax in (("x", gx, xc, (0, 1)), ("y", gy, yc, (0, 2)),
                           ("z", gz, zc, (1, 2))):
        spread = np.abs(g - ref.reshape([-1 if i == 2 else 1 for i in range(3)]
                                        if nm == "x" else
                                        [1, -1, 1] if nm == "y" else [-1, 1, 1])).max()
        step = np.diff(ref).min()
        if spread > 0.01 * step:
            sys.exit("mesh is not axis-separable in %s: spread %.3g vs step %.3g"
                     % (nm, spread, step))
    coords = (zc, yc, xc)
    print("mesh     %d x %d x %d = %d cells (x=scan, y=vertical, z=transverse)"
          % (nx, ny, nz, n))
    vol = read_field(os.path.join(ctl, "0", "V"), 1)
    if vol is None:
        sys.exit("no V field -- run: postProcess -func writeCellVolumes -time 0")

    ts = times(ctl)
    ts = ts[:: max(1, len(ts) // a.max_snaps)][:a.max_snaps]
    print("family   %s" % fam)
    print("control  %s" % a.control)
    print("snapshots %d of %d, transient dropped: t = %s .. %s\n"
          % (len(ts), len(times(ctl)), ts[0], ts[-1]))

    res, dist = {}, {}
    for c in cases:
        cd = os.path.join(fam, c)
        rows = [r for r in (snapshot(cd, t, coords, vol) for t in ts) if r]
        if not rows:
            print("  %-14s no usable snapshots" % c)
            continue
        res[c] = rows
        if c != a.control:
            dist[c] = [r for r in (compare(cd, ctl, t, coords, vol) for t in ts) if r]
        print("  %-14s %3d snapshots, median liquid %.4f mm3"
              % (c, len(rows), med_iqr(rows, "V_liq_mm3")[0]))

    hdr = ("case", "nsnap", "speed", "omega", "fx", "fy", "fz", "up_frac", "Qpos")
    print("\n%-14s %5s %9s %9s %6s %6s %6s %8s %6s" % hdr)
    print("%-14s %5s %9s %9s %6s %6s %6s %8s %6s"
          % ("", "", "m/s", "1/s", "-", "-", "-", "-", "-"))
    for c in [a.control] + [x for x in res if x != a.control]:
        if c not in res:
            continue
        r = res[c]
        print("%-14s %5d %9.4f %9.3g %6.3f %6.3f %6.3f %8.3f %6.3f" % (
            c, len(r), med_iqr(r, "speed")[0], med_iqr(r, "omega")[0],
            med_iqr(r, "fx")[0], med_iqr(r, "fy")[0], med_iqr(r, "fz")[0],
            med_iqr(r, "up_frac")[0], med_iqr(r, "Qpos")[0]))

    if a.control in res:
        b = res[a.control]
        keys = ("fy", "fz", "omega", "speed", "up_frac", "up_flux", "Qpos")
        print("\npaired vs %s -- ratio (Wilcoxon signed-rank p on time-matched snapshots)"
              % a.control)
        print("%-14s" % "case" + "".join("%16s" % k for k in keys)
              + "%8s %8s" % ("D_cos", "D_L2"))
        # Bonferroni over the 7 metrics x however many cases are compared
        ncmp = max(1, (len(res) - 1) * len(keys))
        for c in sorted(x for x in res if x != a.control):
            cells = ""
            for k in keys:
                ratio, pv, n = paired(res[c], b, k)
                if not np.isfinite(ratio):
                    cells += "%16s" % "-"
                    continue
                star = ("**" if pv < 0.05 / ncmp else ("*" if pv < 0.05 else ""))
                cells += "%16s" % ("%.2fx p%.3f%s" % (ratio, pv, star))
            dc = med_iqr(dist[c], "D_cos")[0] if dist.get(c) else np.nan
            dl = med_iqr(dist[c], "D_L2")[0] if dist.get(c) else np.nan
            print("%-14s%s%8.3f %8.3f" % (c, cells, dc, dl))
        print("  ** survives Bonferroni for %d comparisons (alpha=%.5f); * nominal only."
              % (ncmp, 0.05 / ncmp))

    if a.chaos:
        print("\nchaos bound: 0.1 T vs 0.2 T, same plane (a field change of the same kind)")
        print("%-20s %8s %8s" % ("pair", "D_cos", "D_L2"))
        for pl in ("xy", "xz", "yz"):
            lo = os.path.join(fam, "rmf-%s-0.1T" % pl)
            hi = os.path.join(fam, "rmf-%s-0.2T" % pl)
            if not (os.path.isdir(lo) and os.path.isdir(hi)):
                continue
            rows = [r for r in (compare(hi, lo, t, coords, vol) for t in ts) if r]
            if rows:
                print("%-20s %8.3f %8.3f" % ("rmf-%s 0.2 vs 0.1T" % pl,
                                             med_iqr(rows, "D_cos")[0],
                                             med_iqr(rows, "D_L2")[0]))

    print("\nVorticity/Q on the liquid mask eroded by one cell (interface shear is numerical).")
    print("fx+fy+fz = 1: the share of enstrophy about each axis. fz = circulation in the")
    print("scan-vertical plane, the plane a pore must ride to reach the surface.")
    print("up_frac/up_flux describe where the FLOW goes, not where the force points.")


if __name__ == "__main__":
    main()
