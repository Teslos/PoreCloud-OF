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

### 5.5 Parallel verification — 18 ranks

Run on the 0.2 T azimuthal melt-pool case, `scotch`, 18 subdomains, against the identical
serial case. What passed outright:

- **CSV lands in the case root**, not `processor0/` — the `globalPath()` fix of §7.1 holds
  under MPI.
- **No lost or duplicated bubbles.** 181 injected in both runs; `PoreID` unique; the
  gathered list is identically ordered on every rank.
- **The solver fields agree.** `T` and `epsilon1` sampled at the host cell match serial to
  ~1e-4%, so any force difference is the library's, not the flow solution's.

Two real bugs surfaced and were fixed here rather than in production:

1. **The injector max-reduced the cell indices.** `injectorCells_`/`TetFaces_`/`TetPts_`
   are *local* mesh indices and `validInjection()` keys on them being `-1` on non-owning
   ranks. Reducing them made every rank claim every bubble — injecting each one 18 times at
   unrelated cells. `CellZoneInjection` deliberately does not reduce them; neither do we
   now. Velocity and diameter follow the same rule, drawn by the owning rank.
2. **The `maxParcels` cap truncated by prefix.** `setSize(n)` keeps the first `n`, and the
   gathered list is rank-ordered, so the whole quota was filled from the lowest ranks —
   biasing injection toward whatever region those ranks own. Replaced with even striding.

### What does not hold in parallel: sphere-averaged sampling

`fvMesh::cellCells()` does not cross a `processorPolyPatch`, so a bubble straddling a
decomposition plane samples only its local side. 181 bubbles at identical fixed positions
(d/dx 3.1–5.4, ~50 cells per stencil), matched one-to-one to within 1e-11 m:

| | mean coverage | min coverage |
|---|---|---|
| serial | 0.980 | 0.893 |
| 18 ranks | 0.893 | 0.317 |

| parallel coverage | n | mean \|dF\|/F | max |
|---|---|---|---|
| ≥ 0.89 (serial floor) | 130 | 1.79% | 9.14% |
| 0.70 – 0.89 | 24 | 13.05% | 71.76% |
| 0.50 – 0.70 | 21 | 9.45% | 29.90% |
| < 0.50 | 6 | 31.20% | 75.08% |
| **all** | **181** | **5.15%** | **75.08%** |

A control run with centroid sampling — no stencil at all — differs by only 0.54% mean
(15.1% max), which is the `phiE` Poisson solve converging differently under decomposition.
That isolates the remaining error as truncation, not round-off.

**The mitigation, and its limit.** `sphereAverage` now returns the collected volume and
the reporter writes `coverage()` per bubble as a `Coverage` column. In serial nothing fell
below 0.89 — that floor is cell centres failing to tile a sphere exactly, not truncation —
so a parallel bubble materially below it was cut by a processor boundary. Discarding those
keeps 72% of the sample and restores serial-like accuracy: mean 5.15% → 1.79%, max 75% →
9.1%. The floor is resolution dependent; take it from a serial run of the same case rather
than assuming 0.89. `LeenovKolinForce` emits a `WarningInFunction` at construction when
`useSphereAverage` is on in a decomposed run.

This is a reporting fix, not a physics fix. A real one needs remote stencil values, via
`mapDistribute` (as `extendedCentredCellToCellStencil` does) or a layered halo exchange (as
`PoreExtract.C:322-382` does for its CCL). Until then: `useSphereAverage false` is exact in
parallel, and sphere averaging is for serial or coverage-filtered work.

### 5.6 Production run — truncation is a visible force discontinuity

The production-length run that §5.5 called for: 0.2 T azimuthal case, t = 0.001 → 0.002
(the dict's full 1 ms injection window), 18 scotch ranks, 8885 s wall, 8600 steps,
reported every 25 steps. 396 snapshots, 457203 rows. 2892 parcels at the end — 2315 still
active in liquid, 577 entrapped as pores by the solidification front (20%). No duplicated
parcel ids in any snapshot, so the injector fix of §5.5 holds at scale.

Coverage over all 457203 samples: mean 0.877, median 0.962, min 0.098, with **34.1% below
the 0.89 serial floor**. That independently reproduces the 38% truncation fraction §5.5
measured from 181 fixed bubbles, now from a freely-moving sample 2500× larger — so the
one-third figure is a property of the decomposition, not of that particular bubble set.

**The open question from §5.5 is answered: yes, truncation is visible along a trajectory.**
Comparing consecutive samples of the same active bubble, the relative change in |F| is
4.1× larger when its coverage jumps (>0.10) than when coverage is steady (<0.01). That
raw ratio is confounded — a bubble whose stencil changed also moved further, so its fields
genuinely changed more. Controlling for the distance actually travelled between samples:

\verbatim
  displacement    steady n   median    jump n   median    ratio
    < 0.05 dx        6168     2.32%       151   11.12%     4.8x
    0.05-0.2 dx     40716     3.97%      4742   10.68%     2.7x
    0.2-0.5 dx      42845     8.24%     16875   14.11%     1.7x
    0.5-1 dx        20905    18.14%     19093   28.44%     1.6x
    > 1 dx          12466    45.47%     28148   65.37%     1.4x
\endverbatim

The effect survives the control and is strongest exactly where it should be. A bubble that
moved less than 1/20 of a cell should see almost no change in force — and it does not, when
its stencil changes: 11.1% against 2.3%. That is the artefact, cleanly separated from
physical field variation. The ratio decays as displacement grows because real field
variation increasingly dominates the same fixed stencil error.

So the discontinuity is real and sharp, not a slow bias: a bubble crossing a decomposition
plane takes a step change in reported force. For trajectory work this matters more than the
7% mean of §5.5 suggests, because the error is not smooth in time.

**Still not checked:** serial-vs-parallel trajectory agreement at production length. A
serial run of this window is ~18× the wall time, so it needs to be worth the machine.


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
- **Parallel verification** (plan test 6) is done — §5.5 for injection, CSV path and force
  accuracy, §5.6 for the production-length run. Two bugs found and fixed. Truncation is
  confirmed as a step change in force at decomposition planes, not a smooth bias. The one
  remaining piece is serial-vs-parallel trajectory agreement at production length, which
  costs ~18× the wall time of the parallel run.
- **Remote stencil values for `sphereAverage`** (§5.5). Without them, sphere-averaged
  sampling in parallel is accurate only for bubbles above the serial coverage floor. This
  is the one known correctness gap in the library.
- Two-way coupling, via the `fvOptions(rho, U)` hook already present in `UEqn.H`.
