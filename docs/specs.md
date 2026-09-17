# PoreCloud-OF — design notes and validation record

## 1. Why this exists

`AZIMUTHAL_STUDY_SUMMARY.md` §5 reduces pore expulsion to a drift budget:

> escape requires `v_drift · t_liquid > pore depth (~200 µm)`,
> with `v_drift = 2·f·a²/(9µ)`

Every conclusion in that study rests on this one analytic line, applied to a handful of
hand-placed pores of a single radius. Its own stated caveats — *"single realization
(n = 1); pore-count differences are within stochastic noise"* — are exactly what a
Lagrangian cloud removes. And `v_drift ∝ a²` is the strongest lever in the whole study,
yet the resolved-VoF workflow can only test one radius per run.

This library replaces the hand calculation with an integrated trajectory, for a
population, inside the solver, using the same field sampling PoreTracker already uses.

## 2. Architecture

One-way coupled: the cloud only *reads* solver fields. That is what allows it to run as a
`functionObject` against an unmodified `laserbeamFoam` binary, rather than as a solver
fork.

```
laserbeamFoam (unmodified)
  └ Time::run()  ──►  functionObjects::execute()
                        └ poreCloudFunctionObject::execute()
                            ├ mu = rho*nu            (mu is NOT registered by the solver)
                            └ basicKinematicCloud::evolve()
                                ├ forces_.cacheFields(true)
                                │   ├ leenovKolin      → J_MHD (registry) + B (rebuilt)
                                │   └ thermocapillary  → gradT (registry)
                                ├ injectors_.inject()  → meltPoolInjection
                                ├ CloudType::move()    → tracking, parallel transfer
                                │   └ postMove         → solidificationCapture
                                └ functions_.postEvolve → poreForceReport (CSV)
```

### Timing

`laserbeamFoam.C` contains no explicit `functionObjects().execute()`; the hook is inside
`Time::run()`, which is the `while` condition of its time loop. So `execute()` fires once
per timestep **at the top of the loop, before `++runTime`**. `runTime.value()` is the time
of the just-completed step, and every field holds converged end-of-PIMPLE values for that
step — a self-consistent set to advect bubbles through.

Consequence on a **restart**: `J_MHD` is `NO_READ`, so on the first pass after a restart it
is still zero (it is filled during that step's `UEqn.H`). The first step therefore sees no
EM force. Subsequent steps are correct. Immaterial for runs of more than a few steps.

Sub-timestep coupling is not available to a functionObject. The cloud sub-cycles
internally against `maxCo`, which is the right granularity here.

### What had to be rebuilt rather than looked up

`B` is **not registered**. `laserbeamFoam` builds `Bfield_MHD` as a function-local
`volVectorField` inside `UEqn.H`, so it is destroyed before any functionObject runs —
unlike `J_MHD`, `LorentzForce` and `exclusionForce`, which are registered and persist.
`poreCloudMagneticField` therefore re-reads `transportProperties/MHD` and rebuilds `B`
with the solver's own expressions (both `azimuthal` and `uniform` DC/RMF branches,
including the read-time normalisation of `staticField`).

`mu` is **not registered** either: `incompressibleTwoPhaseMixture::mu()` returns a
`NO_REGISTER` temporary. It is formed as `rho*nu` and refreshed every step, because `nu` is
re-blended by the solver's Arrhenius model on every PIMPLE iteration.

## 3. Finite-size sampling

Melt-pool cells are 6–9 µm (`blockMeshDict` fine zones: ~7.1 µm x, ~6.25 µm y, ~9.4 µm z).
The pores of interest are 10–70 µm across — up to ~10 cells. The point-particle assumption
is violated by construction.

`poreCloudSphereAverage.H` breadth-first walks `mesh.cellCells()` from the bubble's host
cell and returns the **volume-weighted mean** of a field over cells within the bubble
radius. The force model multiplies by the bubble's analytic volume `V = 4/3πa³`, so the
product equals PoreTracker's integral `F = -k·Σ(V_c·f_c)`.

A *mean* rather than a raw sum is returned deliberately:

- it degrades gracefully across processor boundaries (`cellCells()` does not cross a
  `processorPolyPatch`, so the walk truncates; normalising keeps the result a correct
  average over the visible part instead of silently losing force);
- it reduces exactly to the host-cell value when the radius is sub-cell, so one code path
  serves both regimes — which is also how `useSphereAverage false` is implemented
  (radius 0).

`dOverDx` is written per bubble per step so the regime is visible in the data. Measured in
the melt-pool smoke test: 1.17–8.15, mean 4.60.

**This is an approximation, not a derivation.** It is defensible and consistent with the
Eulerian tool, but test 3 below is what would bound its error empirically.

## 4. The exclusion-force coefficient

Two implementations in `mhd-laserbeamfoam-solver` disagree by a factor of two:

| Branch | Formula | Location |
|---|---|---|
| `azimuthal-field` (live) | `-1.5·(1-α₁)·(J₀×B)`, `J₀ = fvc::average(fvc::interpolate(J_MHD))` | `UEqn.H:185-202` |
| `fixes-solver` (Chang) | `-0.75·smoothed(LorentzForce)·(1-α₁)·nearLiquid` | `UEqn.H:61-76` |

Commit `0340af85` kept `-3/2`, rejecting `-3/4` as not matching Leenov–Kolin; Chang cites
Takahashi & Taniguchi (2003) for `-3/4`. Both appear in the inclusion-removal literature.

This library does not adjudicate. `exclusionCoeff` defaults to `-1.5` to match the live
Eulerian solver and is a dictionary entry so both can be run. **It scales `v_drift`
linearly and therefore scales every expulsion conclusion**; worth settling from the
primary sources independently of this work.

### No double counting

The Eulerian `exclusionForce` is gated by `(1-α₁)`. Lagrangian bubbles never modify
`alpha.metal`, so `α₁ = 1` at a bubble and that term is identically zero there.

Better: `J_MHD` sampled at a bubble **is** the undisturbed `J₀`, because no void has been
carved out of the conductivity field. The Eulerian version has to reconstruct `J₀` by
smoothing `J` out of the void it does carve, which underestimates for deep voids.

## 5. Validation record

All runs on OpenFOAM v2506, `laserbeamFoam` from branch `azimuthal-field` (binary
unchanged, Aug 21 08:23).

### 5.1 Terminal velocity — the acceptance test

`cases/driftValidation`: 200 µm cube, 25³ uniform cells, all metal, T = 1000 K uniform,
laser off, MHD off, `gravity + sphereDrag + virtualMass`. 200 bubbles at d = 35 µm.
Terminal velocity is analytic.

```
mu = rho*nu = 2670 * 6.71677e-8*exp(2078.94/1000) = 1.4340e-3 Pa.s
Re = 0.081                       (<< 1, Stokes regime)

v_Stokes            = 1.2426 mm/s
v_Schiller-Naumann  = 1.2103 mm/s       <- what sphereDrag implements
measured            = 1.2058 mm/s       -0.38 % vs SN, -2.97 % vs Stokes
sd over 198 bubbles = 4.8e-10 mm/s
mean Ux, Uz         = 1.2e-18, 2.5e-19 mm/s    (machine zero)
```

Approach to terminal is monotone, reaching 99 % at ~320 µs — consistent with the
added-mass-dominated response time `τ = (ρ_p + ½ρ_c)d²/(18µ) ≈ 63 µs`.

This validates buoyancy, drag, added mass, the velocity integration, and the whole
functionObject → cloud → evolve path. It is the same balance as
`v_drift = 2fa²/(9µ)` with `f = (ρ_m - ρ_g)g`.

### 5.2 EM force — orthogonality

`F = C·V·(J₀ × B)` must be exactly perpendicular to `B`. Tested on the melt-pool case with
a **uniform** DC field (B constant everywhere, so no spatial-variation confound):

```
rows with non-zero F_em = 1779
max |cos(F_em, B)|      = 0.000e+00
max |Fz|/|F|            = 0.000e+00     (B along +z)
```

Exact to machine precision. Validates the cross product, the coefficient path and the B
reconstruction.

With the **azimuthal** field the same test gives `max|cos| = 3.8e-2`. That is not an
error: `B` varies steeply across a 35 µm bubble near the axis (`1/√(r²+r_c²)`), so the
force is perpendicular to `⟨B⟩` over the sphere, not to `B` at the centroid. Sampling at
the centroid instead gives `max|cos| = 0.64`, because the *cell-centre* `B` used by the
code then differs in direction from the `B` at the bubble position used by the check.

### 5.3 Force magnitude vs the published figure

Melt-pool case, uniform 0.2 T, EM force density `f = |F_em|/V` at bubbles:

```
median f    = 2.645e5 N/m3
study cites = 2.9e5   N/m3      (AZIMUTHAL_STUDY_SUMMARY.md section 5)
```

Independent agreement to ~9 % on the force magnitude.

### 5.4 A discrepancy in the quoted drift velocity

The force density agrees, but the **drift velocity quoted in §5 of the study does not
follow from its own formula and its own numbers**. With `f = 2.9e5 N/m³`, `a = 35 µm`:

| T [K] | µ = ρν [Pa·s] | `v_drift = 2fa²/9µ` [mm/s] |
|---|---|---|
| 900 | 1.807e-3 | 10.9 |
| 1000 | 1.434e-3 | 13.8 |
| 1500 | 7.171e-4 | 27.5 |
| 2000 | 5.071e-4 | 38.9 |
| 2743 (T_vap) | 3.827e-4 | 51.6 |
| — | 3.289e-4 | **60 (as quoted)** |

The quoted 60 mm/s implies µ = 3.29e-4 Pa·s, i.e. ν = 1.23e-7 m²/s — **below the solver's
own Arrhenius viscosity at any temperature, including the vaporisation point**. At
melt-pool temperatures away from the keyhole (900–1500 K) the formula gives 11–28 mm/s,
2–5× lower.

This does not overturn the study's conclusion — it concluded expulsion *fails* at 0.2 T,
and a lower drift strengthens that. But it materially changes the §6 extrapolation of what
field strength or dwell time would be needed: the gap is nearer 10–20× than the stated
3–4×. Integrating the trajectory rather than estimating it is precisely what removes this
class of uncertainty.

*(Stated as an arithmetic observation. The viscosity actually used for the quoted figure
is not recorded in the summary, so this is a reconciliation, not an accusation.)*

## 6. Deliberate omissions

**Drag / added mass / lift are not broken out in the CSV.** They are evaluated inside
`KinematicParcel::calcVelocity` from `trackingData` that is only valid mid-move and is
stale by `postEvolve`. Recomputing them in the reporter would duplicate each drag
correlation and risk silent divergence from what was actually integrated. `Uslip` is
reported instead, from which any drag law can be evaluated in post-processing against the
same slip velocity the solver used.

**All body forces return `Su`, never `Sp`.** In `KinematicParcel.C:207-229` the "proper
splitting" that would honour an implicit non-coupled force is commented out —
`Fncp.Sp()` is silently discarded.

**Capture is irreversible by default.** Releasing a captured bubble would mean re-deciding
where the front is relative to a body that is no longer being advected. `allowRemelt true`
enables it.

## 7. Two PoreTracker bugs not reproduced

1. **CSV path under MPI.** `PoreExtract` writes to `mesh_.time().path()`, which resolves to
   the *processor-local* case directory — so its CSV lands in `processor0/`, not the case
   root. This library uses `time().globalPath()`.

2. **Write cadence.** OpenFOAM only wraps a functionObject in `timeControl` when
   `writeControl`/`executeControl` are present in its sub-dict
   (`functionObjectList.C:1066,1126`). PoreTracker's shipped dict sets neither, so despite
   its comments it writes a row per pore *per timestep*. `poreForceReport` implements its
   own explicit `writeControl` (`writeTime` | `timeStep`), and the shipped dicts declare it.

Also fixed: the CSV is opened lazily on first write and in **append** mode, so a restart
continues the series rather than truncating it, and a cloned function object does not end
up holding a null stream.

## 8. Remaining work

- **Test 3 from the plan is not yet done**: cross-check against the resolved-VoF results by
  injecting bubbles at the `pore_coords.json` sites and comparing trajectories with the
  Jul 14 `poreTracker1_pores.csv`. Those results predate the Aug 21 merge (`cd729daf`), so
  that case needs re-running on the current binary first to get a clean baseline. This is
  the test that would bound the finite-size approximation empirically.
- **Parallel verification** (plan test 6) has not been run: serial-vs-parallel trajectory
  agreement, CSV landing in the case root, no lost/duplicated bubbles across processor
  boundaries, and whether `sphereAverage` truncation at decomposition planes is visible as
  a discontinuity in force along a trajectory. The code paths are written for MPI
  (`globalIndex` gather in the injector, `Pstream::gatherList` in the reporter) but are
  untested under it.
- Two-way coupling, via the `fvOptions(rho, U)` hook already present in `UEqn.H`.
