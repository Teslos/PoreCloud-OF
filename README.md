# PoreCloud-OF

Lagrangian pore/bubble cloud for MHD melt-pool simulations, for OpenFOAM v2506 (ESI).

Injects gas bubbles into the liquid region of a laser melt pool and integrates the force
balance acting on them — buoyancy, drag, added mass, fluid acceleration, thermocapillary
migration, and the electromagnetic Leenov–Kolin exclusion force — to quantify how an
applied magnetic field moves, traps or expels porosity during laser welding.

It attaches as a runtime-loaded `functionObject` to an **unmodified `laserbeamFoam`
binary**. No solver fork, no recompile of the solver.

## Relationship to PoreTracker-OF

The two are complements, not alternatives.

| | PoreTracker-OF | PoreCloud-OF |
|---|---|---|
| Paradigm | Eulerian, observational | Lagrangian, predictive |
| Finds pores by | connected-component labelling on `alpha.metal` | injecting them |
| Sees | voids the VoF scheme already resolves | bubbles of any size, resolved or not |
| Forces | samples `LorentzForce` at pore cells | integrates a full force balance |
| Population | whatever the run happened to produce | a controlled size distribution |

`PoreCloud-OF` writes a CSV whose first fifteen columns are exactly PoreTracker's, so the
existing analysis chain (`PoreStats.py`, `visualize_poretracker.py`,
`plot_pore_trajectories_csv.py`, `compare_poretracker.py`) reads its output unmodified.

It replaces the older manual workflow in `mhd-openfoam` — `findPoreSites.py` →
`injectPore` → `setFields` → restart — which could only place a handful of hand-positioned
pores of a single radius, all at one instant, and required the mesh to resolve them.

## Requirements

- OpenFOAM v2506 (ESI). Uses `liblagrangianIntermediate`'s cloud, force, injection and
  cloud-function-object machinery.
- A solver that registers `U`, `rho`, `nu`, `T`, `epsilon1`, `alpha.metal`, `gradT` and
  (for the EM force) `J_MHD` — i.e. `laserbeamFoam` from the `mhd-laserbeamfoam-solver`
  repo, branch `azimuthal-field` or later.

## Build

```bash
source /usr/lib/openfoam/openfoam2506/etc/bashrc
./Allwmake          # -> $FOAM_USER_LIBBIN/libPoreCloud.so, $FOAM_USER_APPBIN/poreCloudReplay
./Allwclean
```

## Use

Add to `system/controlDict`:

```cpp
functions
{
    poreCloud
    {
        type            poreCloud;
        libs            ( "libPoreCloud.so" );
        executeControl  timeStep;
        executeInterval 1;
    }
}
```

and copy `cases/meltPool-azimuthal/constant/poreCloudProperties` into your case's
`constant/`. Sub-models are read from `constant/<cloudName>Properties`
(`poreCloudProperties` by default).

Output: `<caseDir>/poreCloud_pores.csv`, plus `<time>/lagrangian/poreCloud/` for ParaView.

## Replaying a finished run (`poreCloudReplay`)

The cloud is one-way coupled — it reads the melt fields, the melt never sees the
parcels — and measured cost is ~0.1% of runtime against ~100% for the CFD. Re-solving a
2.5-hour `laserbeamFoam` run just to try a different cloud parameter (`exclusionCoeff`,
injection rate, size distribution, ...) is therefore pure waste.

`applications/poreCloudReplay` is a standalone application that drives the same,
unmodified `poreCloud` functionObject against the already-stored fields of a finished
run instead. It reads a case's own mesh and controlDict as any OpenFOAM application
would, discovers the stored snapshot times, then advances with its own (much finer)
timestep, linearly interpolating `U`, `J_MHD`, `epsilon1`, `T` and `alpha.metal` in time
between the two stored snapshots bracketing the replay clock. Disk is touched only when
the replay clock crosses into a new bracket — reading a full timestep of fields on every
sub-step would cost as much as the CFD it replaces. `rho` and `nu` are not stored on
disk (`laserbeamFoam` never writes them) and are re-derived every step from the
interpolated `alpha.metal`, using the same formulas `createFields.H` /
`laserbeamFoam.C` use (mass-weighted Newtonian blend, or the Arrhenius `nu(T)` override
when the case's `transportProperties` defines it).

```bash
cd <finished-case>      # already has stored time directories and
                         # constant/poreCloudProperties + the functionObject entry
poreCloudReplay          # serial
mpirun -np N poreCloudReplay -parallel   # or against a decomposed case
```

Set `startTime`/`endTime`/`deltaT`/`adjustTimeStep` in `system/controlDict` as for any
solver; `startTime` must fall inside the stored snapshot range (replay cannot
extrapolate beyond it). Serial is preferred when available, for the same
`cellCells()`-does-not-cross-a-processor-boundary reason noted below for sphere
averaging.

**Limitation to validate before trusting the output**: storage cadence is typically
1e-5 s against a CFD timestep around 1e-7 s — roughly 87x coarser. Linear interpolation
between snapshots that far apart smooths out anything that varies faster than the
storage cadence (fast MHD/turbulent transients, a rotating field with a period
comparable to the cadence, non-monotonic behaviour between two snapshots). Treat replay
output as a fast first look, not ground truth, until checked against a short
directly-coupled run over the same window. See the header comment in
`poreCloudReplay.C` for the full reasoning.

## What this library provides

Most of the machinery is stock OpenFOAM. New here:

| Component | Type | Purpose |
|---|---|---|
| `poreCloud` | functionObject | owns and evolves the cloud against an unmodified solver |
| `leenovKolin` | particle force | EM exclusion force, `F = C·V·(J₀ × B)` |
| `thermocapillary` | particle force | Young–Goldstein–Block Marangoni migration |
| `meltPoolInjection` | injection model | continuous nucleation where the metal is molten |
| `solidificationCapture` | cloud function | freezes a bubble once the front passes it |
| `poreForceReport` | cloud function | PoreTracker-compatible per-bubble CSV |

Reused unchanged: `gravity` (which already includes buoyancy), `sphereDrag` /
`TomiyamaDrag`, `virtualMass`, `pressureGradient`, `TomiyamaLift`, `interface`, and
`dataCloud` / `vtkCloud` for ParaView output.

## Two things to know before trusting the output

**1. The exclusion-force coefficient is contested, and it is a dict entry for that reason.**

`mhd-laserbeamfoam-solver` contains two independently-derived implementations that
disagree by a factor of two — `-3/2·(J₀×B)` on `azimuthal-field`, `-3/4·F_Lorentz` on
`fixes-solver` (Chang, citing Takahashi 2003). Commit `0340af85` kept `-3/2`. Both
coefficients appear in the inclusion-removal literature. Since the coefficient scales the
drift velocity linearly, it scales every expulsion conclusion. `exclusionCoeff` defaults
to `-1.5` to match the live Eulerian solver; set it to `-0.75` to test the alternative.

**2. These bubbles are larger than the cells, deliberately.**

Melt-pool cells are 6–9 µm; the pores of interest are 10–70 µm across. The
point-particle assumption behind Lagrangian tracking does not hold. The EM and
thermocapillary forces therefore volume-average the carrier fields over the bubble sphere
(`poreCloudSphereAverage.H`), reproducing the integral PoreTracker performs over a
resolved pore. The `dOverDx` CSV column reports `d/Δx` per bubble per step so the regime
stays visible — in the smoke tests it runs 1.2–8.2, mean 4.6.

## Validation

`cases/driftValidation` is a quiescent, uniformly molten box with the laser and MHD off,
where the terminal velocity is analytic. Run it with `./Allrun`:

```
v_Stokes             1.2426 mm/s
v_Schiller-Naumann   1.2103 mm/s     <- what sphereDrag should give
measured             1.2058 mm/s     (-0.38 %, 198 bubbles, sd 4.8e-10)
lateral drift        ~1e-18 mm/s     (machine zero, as it must be)
```

See `docs/specs.md` for the full validation record, including the EM-force orthogonality
test and a note on the drift-velocity figure quoted in `AZIMUTHAL_STUDY_SUMMARY.md`.

## Status

One-way coupled: bubbles read the melt and feel every force, but exert no reaction on `U`.
Two-way coupling is a later phase and would use the `fvOptions(rho, U)` hook already
present in `UEqn.H`.

## Author

Toni Ivas
