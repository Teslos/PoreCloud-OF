# Can a magnetic field push pores out of a laser melt pool?

**Short answer: no field configuration tested does.** The best are directionally
neutral; the azimuthal field, which the campaign was built around, is the worst of
those tested and pushes pores *away* from the surface.

This page summarises the evidence, the two methods used, and a methodological result
that affects how any future pore study should be run. It is written to stand alone —
no prior context assumed.

---

## 1. The question

Gas pores form in the melt pool during laser welding and LPBF, and become defects when
they freeze in. If an electromagnetic body force could drive them to the free surface
before solidification, they would escape.

The relevant force is the **Leenov–Kolin exclusion force** on a non-conducting sphere
in a conducting liquid carrying current density **J** in field **B**:

```
F_exclusion = -1.5 · V_bubble · (J × B)
```

Note the sign: it is *antiparallel* to the Lorentz body force **J × B** that acts on the
melt itself. Expulsion means reaching the free surface, which in this geometry is the
**+y** direction.

## 2. Two independent methods

| | what it measures | cost |
|---|---|---|
| **Lagrangian cloud** | force on individual tracked bubbles, integrated along trajectories | ~2.5 h CFD per case |
| **Eulerian survey** | direction of the force field throughout the liquid, per stored snapshot | seconds |

The Eulerian method works because `F_exclusion` is antiparallel to the `LorentzForce`
field the solver already writes — confirmed empirically at `cos(angle) = −0.9993`. So
the direction question needs no particle tracking and no time integration.

## 3. Result

Five configurations run with the Lagrangian cloud, same t = 1→2 ms window, per-case
calibrated injection, 18 MPI ranks:

| configuration | bubbles | **% pushed up** | 95% CI | % that rose | % reaching surface |
|---|---|---|---|---|---|
| azimuthal 0.2 T | 330 | **36.1%** | 30.9–41.2 | 43.6% | 7.9% |
| azimuthal **1 T** | 85 | 32.9% | 22.9–42.9 | **25.9%** | 3.5% |
| DC, B∥x, 0.2 T | 236 | **48.3%** | 41.9–54.7 | 44.9% | 8.9% |
| RMF xz-plane, 0.2 T | 386 | **47.4%** | 42.4–52.4 | 45.6% | 8.8% |
| RMF yz-plane, 0.2 T | 327 | **48.9%** | 43.5–54.3 | 42.8% | **11.0%** |

**Nothing exceeds 50%.** DC-x, RMF-xz and RMF-yz are statistically indistinguishable
from one another and all sit just below neutrality. Azimuthal is the worst
(+12.2 points to DC-x, z = 2.9).

**More field makes it worse, not better.** At 1 T the direction statistic is unchanged
against 0.2 T (z = −0.5) — as expected, since the current here is Seebeck-dominated and
therefore independent of B, so `F = J × B` grows linearly with B but does not change
direction. Yet the fraction of bubbles that rose falls **43.6% → 25.9%** (z = −3.2). A
stronger force pointing the wrong way is worse than a weaker one.

**The gain from abandoning azimuthal is real but modest:** 8.8–11.0% of bubbles reach
the surface against 7.9%, achieved by *removing a downward force* rather than adding an
upward one.

Across all 35 archived cases the Eulerian survey adds two structural findings:

- **Orientation dominates strength.** 0.1 T vs 0.2 T moves RMF results by under 1.5
  points in every matched pair, and 0.2 T vs 1 T does not move the direction statistic
  at all. The current is Seebeck-dominated, hence B-independent, so `F = J × B` scales
  linearly with B while its direction is set by the thermoelectric current structure.
- **B∥y can never help.** Since `J × B ⊥ B`, a field along y gives `F_y ≡ 0` exactly —
  verified as zero in all 57,306 cells carrying force. It is not a poor configuration;
  it is structurally incapable of vertical force.

## 4. Two traps for anyone repeating this

### 4.1 Injected resolved pores do not survive

A 35 µm-radius gas void injected into the melt collapses or is absorbed **within ~50 µs**:

```
t = 500 us (injection)   4 isolated voids   [437, 437, 435, 435] cells
t = 550 us               7 isolated voids   [12, 10, 6, 2, 1, 1] cells
```

Nothing holds such a void open against ambient pressure in this solver. Anything
measured on "injected pores" more than a few tens of microseconds after injection is
measuring something else — usually the keyhole wall, since that is the other large gas
region adjacent to melt.

**Always verify survival** with a connected-component check on `alpha.metal` a few write
intervals after injection. It is a few lines and it changes conclusions.

A related timing trap: at t = 100 µs the pool holds **105 liquid cells** while a pore
occupies **398**. Early in a run there is physically nowhere to place a pore that is not
the keyhole, regardless of the coordinates requested.

### 4.2 A volume-averaged force field is not what bubbles experience

The Eulerian survey is roughly unbiased but imprecise. Against four tracked-bubble runs
its deltas were **−1.1, +8.6, +4.9 and −5.0 points** — mean +1.8, standard deviation 6.0,
and the signs differ. Bubbles are born in liquid, advected by the flow and frozen where
they solidify, so they do not sample the pool the way a volume average assumes, and the
difference does not fall the same way in every case.

A spread of 6 points cannot separate configurations differing by 1–2, which is what the
three best ones do. Use the survey for the coarse split — azimuthal against the rest,
which it gets right by a wide margin — and the cloud for a number.

## 5. Reconciling the earlier azimuthal study

An earlier study reported an up-fraction of **0.78** for the azimuthal field — the
opposite conclusion. Both results are correct; they measure different things.

That study's mask, `gas & dilate(melt,2) & T>870`, does not exclude the keyhole.
Decomposed over its own analysis window, it contains **6076 keyhole-connected cells
against 151 isolated-pore cells — 97.6% keyhole wall**. Its pores had long since
collapsed (§4.1), leaving the keyhole as the only gas adjacent to melt.

The azimuthal field is *defined* to circulate around the keyhole axis, so a coherent
surface-ward force on that wall is close to a restatement of the field's construction.

**The 0.78 is valid as a keyhole-wall measurement and should be labelled as one.** It is
not evidence that the field expels pores.

## 6. What would change the answer

- **A stronger or differently-shaped field.** Orientation matters more than strength, and
  nothing tested reaches 50%, so the search space is not exhausted.
- **Two-way coupling.** The cloud is one-way coupled: bubbles read the melt, the melt
  never sees them. Valid at the measured void fractions (~1–3%), but it means bubbles
  cannot collectively alter the flow that carries them.
- **Where pores actually nucleate.** Bubbles here are seeded throughout the liquid. If
  formation is strongly localised — at the keyhole, say — the relevant force is the local
  one, and the population statistic would need reweighting.

## 7. Reproducing

```bash
source /usr/lib/openfoam/openfoam2506/etc/bashrc
# analysis scripts need the project venv (system python3 has no pandas)
VENV=~/.venvs/mhd-paper-analysis/bin/python3

# Eulerian survey across cases - seconds
$VENV scripts/poreCloudFieldSurvey.py <case> [<case> ...] -o summary.csv

# Full expulsion budget from a finished cloud run
$VENV scripts/poreCloudReport.py <case>/poreCloud_pores.csv -o report.pdf

# Compare several finished runs
$VENV scripts/poreCloudCompare.py label=<case> ... -o compare.pdf
```

Injection **must be calibrated per case** — `rate` is a nucleation rate per unit *liquid*
volume per second, so the resulting count is `rate · V_melt · T`. An uncalibrated value
once produced a bubble population occupying **473% of the melt volume**, which nothing in
the solver objects to because the coupling is one-way. `PoreForceReport` now reports void
fraction every write and warns past 5%; the resolved-VoF runs sit at 0.8–2.5%.

## 8. Where things live

| | |
|---|---|
| library, analysis scripts, full technical record | `PoreCloud-OF` (`docs/specs.md`) |
| per-case injection calibration | `~/results/poreCloud-calibration/` |
| 35-case Eulerian survey | `~/results/poreCloud-survey/` |
| Lagrangian runs | `~/results/poreCloud-{0.2T-azimuthal-18rank,dc-bx-0.2T,rmf-xz-0.2T}/` |
| pore-survival experiments | `~/results/azimuthal-pores-{rerun,offaxis}/` |

`docs/specs.md` §5.11 and §5.12 carry the detailed numbers, uncertainties and the
reasoning behind each conclusion.
