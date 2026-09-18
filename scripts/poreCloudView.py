#!/usr/bin/env python3
"""Build the standard poreCloud view in ParaView: pores over the melt pool.

Two ways to use it.

  Regenerate the image and a loadable state file, headless:

      pvpython poreCloudView.py --case <case>/<name>.foam \\
                               --pores <dir>/pores.pvd \\
                               --png out.png --state out.pvsm

  Or build the same pipeline inside a running ParaView GUI - open
  Tools > Python Shell (the shell runs on the server you are connected to) and:

      exec(open('/path/to/poreCloudView.py').read().split('if __name__')[0])
      build(case='<case>.foam', pores='<dir>/pores.pvd')

The view: the metal/gas interface and the liquid contour, both nearly
transparent so the interior is visible, with bubbles as spheres sized by their
real diameter and coloured by the vertical component of the Leenov-Kolin
force.  Blue is downward - the direction that works against expulsion - so the
picture reads as the specs 5.8 result rather than as decoration.
"""

import argparse

from paraview.simple import *  # noqa: F403


# Nearly transparent: bubbles live *inside* the pool, so an opaque interface
# hides every one of them.  This is the single most common reason the cloud
# looks empty in a fresh view.
IFACE_OPACITY = 0.16
LIQUID_OPACITY = 0.22

# Symmetric about zero so the colour map reads as a direction, not a magnitude.
FORCE_RANGE = 2e-8


def build(case, pores, time=None, view=None):
    """Create the pipeline. Returns (view, time_used)."""
    v = view or CreateRenderView()
    v.ViewSize = [1400, 900]
    v.Background = [1, 1, 1]
    v.UseColorPaletteForBackground = 0

    foam = OpenFOAMReader(FileName=case)
    # Without this a parallel run that was never reconstructed shows exactly
    # one timestep, which looks like missing data rather than a wrong setting.
    foam.CaseType = 'Decomposed Case'
    foam.UpdatePipelineInformation()
    foam.MeshRegions = ['internalMesh']
    foam.CellArrays = ['alpha.metal', 'T', 'epsilon1']

    times = list(foam.TimestepValues)
    t = times[-1] if time is None else min(times, key=lambda x: abs(x - time))

    c2p = CellDatatoPointData(Input=foam)

    iface = Contour(Input=c2p, ContourBy=['POINTS', 'alpha.metal'],
                    Isosurfaces=[0.5])
    di = Show(iface, v)
    di.Opacity = IFACE_OPACITY
    ColorBy(di, None)
    di.AmbientColor = di.DiffuseColor = [0.45, 0.50, 0.58]

    liquid = Contour(Input=c2p, ContourBy=['POINTS', 'epsilon1'],
                     Isosurfaces=[0.5])
    dl = Show(liquid, v)
    dl.Opacity = LIQUID_OPACITY
    ColorBy(dl, None)
    dl.AmbientColor = dl.DiffuseColor = [0.85, 0.55, 0.15]

    pvd = PVDReader(FileName=pores)
    calc = Calculator(Input=pvd)
    calc.ResultArrayName = 'F_EM_y'
    calc.Function = 'F_exclusion_Y'

    glyph = Glyph(Input=calc, GlyphType='Sphere')
    glyph.ScaleArray = ['POINTS', 'd_m']   # real bubble size, not a fixed dot
    glyph.ScaleFactor = 0.5
    glyph.GlyphMode = 'All Points'

    dg = Show(glyph, v)
    ColorBy(dg, ('POINTS', 'F_EM_y'))
    lut = GetColorTransferFunction('F_EM_y')
    lut.ApplyPreset('Cool to Warm (Extended)', True)
    lut.RescaleTransferFunction(-FORCE_RANGE, FORCE_RANGE)
    dg.SetScalarBarVisibility(v, True)
    bar = GetScalarBar(lut, v)
    bar.Title = 'EM force, y  (N)'
    bar.ComponentTitle = ''
    bar.TitleColor = bar.LabelColor = [0.1, 0.1, 0.1]

    for src in (foam, pvd, calc, glyph):
        UpdatePipeline(time=t, proxy=src)

    # The animation scene, not just the view, is what a saved state restores.
    # Setting only view.ViewTime gives a .pvsm that loads back at t = 0 - where
    # the cloud is genuinely empty, so the view looks broken on open.
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    scene.AnimationTime = t
    v.ViewTime = t

    v.ResetCamera()
    cam = GetActiveCamera()
    cam.Elevation(-18)
    cam.Azimuth(38)
    v.CameraViewAngle = 28
    Render(v)
    return v, t


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--case', required=True, help='<case>/<name>.foam')
    ap.add_argument('--pores', required=True, help='pore_viz/pores.pvd')
    ap.add_argument('--time', type=float, default=None,
                    help='simulation time (default: last, where the cloud is fullest)')
    ap.add_argument('--png', default=None)
    ap.add_argument('--state', default=None, help='write a .pvsm to load in the GUI')
    args = ap.parse_args()

    v, t = build(args.case, args.pores, args.time)
    print(f'built at t = {t:.6f} s')
    if args.png:
        SaveScreenshot(args.png, v, ImageResolution=[1400, 900])
        print(f'image -> {args.png}')
    if args.state:
        SaveState(args.state)
        print(f'state -> {args.state}')


if __name__ == '__main__':
    main()
