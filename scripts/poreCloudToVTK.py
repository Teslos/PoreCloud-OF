#!/usr/bin/env python3
"""Write poreCloud_pores.csv as a ParaView VTP series.

OpenFOAM already writes the cloud itself into <time>/lagrangian/poreCloud, and
ParaView reads that straight from the decomposed case.  What it does *not*
carry is anything the reporter computed: the Leenov-Kolin force, buoyancy,
stencil Coverage, the active/captured flag.  Those live only in the CSV, and
they are the reason to look at the pores at all.

So this writes one .vtp per reported time plus a .pvd manifest, with every
scalar and vector the CSV carries attached to the points.  Open the .pvd and
ParaView has the timesteps already, at the real simulation times.

    python3 poreCloudToVTK.py poreCloud_pores.csv -o pore_viz

Colour by Active to see capture, by Coverage to see which samples are cut by a
processor boundary, by Fy_N to see which way the EM force is pushing.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

# Point data carried onto every pore.  Vectors first, then scalars.
VECTORS = {
    'F_exclusion': ('Fx_N', 'Fy_N', 'Fz_N'),
    'F_buoyancy': ('Fbuoyx', 'Fbuoyy', 'Fbuoyz'),
    'U': ('Ux', 'Uy', 'Uz'),
    'Uslip': ('Uslipx', 'Uslipy', 'Uslipz'),
}
SCALARS = ['PoreID', 'Volume_m3', 'd_m', 'dOverDx', 'NCells',
           'T_K', 'epsilon1', 'Coverage', 'Active']


def _f(row, key):
    try:
        return float(row[key])
    except (KeyError, ValueError):
        return 0.0


def write_vtp(path, rows, buoy_sign=1):
    """One timestep as a VTK PolyData of vertices."""
    n = len(rows)
    out = ['<?xml version="1.0"?>',
           '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">',
           '  <PolyData>',
           f'    <Piece NumberOfPoints="{n}" NumberOfVerts="{n}">',
           '      <Points>',
           '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">']
    out += [f'          {_f(r, "Cx_m")} {_f(r, "Cy_m")} {_f(r, "Cz_m")}' for r in rows]
    out += ['        </DataArray>', '      </Points>',
            '      <Verts>',
            '        <DataArray type="Int32" Name="connectivity" format="ascii">']
    out.append('          ' + ' '.join(str(i) for i in range(n)))
    out += ['        </DataArray>',
            '        <DataArray type="Int32" Name="offsets" format="ascii">']
    out.append('          ' + ' '.join(str(i + 1) for i in range(n)))
    out += ['        </DataArray>', '      </Verts>', '      <PointData>']

    for name, cols in VECTORS.items():
        sign = buoy_sign if name == 'F_buoyancy' else 1
        out.append(f'        <DataArray type="Float32" Name="{name}" '
                   f'NumberOfComponents="3" format="ascii">')
        out += [f'          {sign*_f(r, cols[0])} {sign*_f(r, cols[1])} '
                f'{sign*_f(r, cols[2])}' for r in rows]
        out.append('        </DataArray>')

    # Magnitudes, so a glyph can be scaled without a Calculator filter first.
    for name, cols in (('F_exclusion_mag', VECTORS['F_exclusion']),
                       ('speed', VECTORS['U'])):
        out.append(f'        <DataArray type="Float32" Name="{name}" format="ascii">')
        out += ['          %g' % (sum(_f(r, c) ** 2 for c in cols) ** 0.5)
                for r in rows]
        out.append('        </DataArray>')

    for s in SCALARS:
        out.append(f'        <DataArray type="Float32" Name="{s}" format="ascii">')
        out += ['          %g' % _f(r, s) for r in rows]
        out.append('        </DataArray>')

    out += ['      </PointData>', '    </Piece>', '  </PolyData>', '</VTKFile>']
    path.write_text('\n'.join(out) + '\n')


def write_series(csv_path, out_dir, buoy_sign=1, stride=1):
    out_dir.mkdir(parents=True, exist_ok=True)
    by_time = defaultdict(list)
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            by_time[float(row['Time'])].append(row)

    times = sorted(by_time)[::stride]
    entries = []
    for i, t in enumerate(times):
        name = f'pores_{i:05d}.vtp'
        write_vtp(out_dir / name, by_time[t], buoy_sign)
        entries.append((t, name))

    pvd = ['<?xml version="1.0"?>',
           '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
           '  <Collection>']
    pvd += [f'    <DataSet timestep="{t}" group="" part="0" file="{n}"/>'
            for t, n in entries]
    pvd += ['  </Collection>', '</VTKFile>']
    manifest = out_dir / 'pores.pvd'
    manifest.write_text('\n'.join(pvd) + '\n')
    return manifest, len(entries), sum(len(by_time[t]) for t in times)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv', type=Path)
    ap.add_argument('-o', '--out-dir', type=Path, default=None)
    ap.add_argument('--stride', type=int, default=1,
                    help='write every Nth reported time (default all)')
    ap.add_argument('--buoyancy-sign', type=int, choices=(-1, 1, 0), default=0,
                    help='0 = auto-detect the pre-2026-09-17 inverted convention')
    args = ap.parse_args()

    sign = args.buoyancy_sign
    if sign == 0:
        import sys
        sys.path.insert(0, str(Path(__file__).parent))
        from poreCloudReport import buoyancy_sign, load_rows
        sign = buoyancy_sign(load_rows(args.csv))
        if sign == -1:
            print('note: CSV predates the buoyancy sign fix; flipping F_buoyancy')

    out = args.out_dir or args.csv.parent / 'pore_viz'
    manifest, nt, npts = write_series(args.csv, out, sign, args.stride)
    print(f'{nt} timesteps, {npts:,} points -> {manifest}')
    print(f'Open {manifest} in ParaView; colour by Active, Coverage or F_exclusion.')


if __name__ == '__main__':
    main()
