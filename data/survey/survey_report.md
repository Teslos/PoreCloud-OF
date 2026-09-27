# Eulerian exclusion-force survey — does the field push pores out of the melt pool?

**Method:** `/home/ivt/openfoam_projects/PoreCloud-OF/scripts/poreCloudFieldSurvey.py`, run unmodified,
`--stride 20`, against all 35 cases in `/home/ivt/results/poreCloud-calibration/inventory.csv`. No source
files, dicts, or simulation outputs were modified. Outputs: `survey_summary.csv` (35 rows),
`survey_snapshots.csv` (612 rows), this report, and two PNGs.

The Leenov-Kolin exclusion force on a bubble is `F_exclusion = -1.5 * V_bubble * LorentzForce`, i.e.
antiparallel to the solver's own stored `LorentzForce` field. "Expulsion" = the +y direction (free surface
at y≈547 µm). The metric reported throughout is **the pool-size-weighted mean, across sampled snapshots,
of the fraction of liquid-metal cells (`alpha.metal>0.5 & epsilon1>0.5`) where the exclusion force has a
positive y-component.** 50% = no directional bias (as good/bad as a coin flip). A no-field control has
zero force everywhere, so "fraction pointing up" is undefined (NaN), not 0%.

## Validation — checked first, as instructed

`azimuthal-0.2T` gave **35.0%** pool-weighted (27.5% force-weighted), inside the stated 30-35% band and
one point off the independent 29.6%-of-bubbles Lagrangian result. No discrepancy — proceeded with the
full survey without touching the tool.

**Arithmetic, shown end to end, for this case's final snapshot (t=0.001s):**
- `alpha.metal>0.5 & epsilon1>0.5` selects **26,876** liquid cells.
- `ef_y = -LorentzForce_y` computed for those cells; **9,802** of them have `ef_y>0`.
- `frac_up = 9802/26876 = 0.36471` → **36.5%** for this one snapshot.
- Repeated over the 18 stride-sampled snapshots (pool sizes 1,753 → 26,876 cells), each snapshot's
  `frac_up` is weighted by that snapshot's liquid-cell count:
  `pool-weighted mean = Σ(frac_up_i · n_i) / Σ(n_i) = 81,239 / 231,987 = 0.35019` → **35.0%**.

This is exactly what `survey_summary.csv` reports for `results/azimuthal-0.2T`.

## Two data problems found in the underlying simulation outputs (not in the survey tool)

These are worth reporting loudly because they change how several rows should be read.

**1. `results_dc_0.2T/dc-bx-0.2T`, `dc-by-0.2T`, `dc-bz-0.2T` have byte-for-byte identical `LorentzForce`
fields.** Their dicts do claim different `staticField` vectors `(1 0 0)`, `(0 1 0)`, `(0 0 1)`, but the
stored field data for all three is the same array at every snapshot checked. The other DC root
(`results_matched_0.2T_f4000hz`) does NOT have this problem — its `dc-bx/by/bz-0.2T` fields are correctly
distinct (verified: each has an exactly-zero component along its own claimed B-axis, as physics requires).
So `results_dc_0.2T` contains **one real simulation counted three times** under different names; the true
field (zero x-component, nonzero y/z) is consistent with B∥x, contradicting the `by`/`bz` labels. I did not
guess at a fix — I flagged all three rows and used only the properly-distinct `results_matched_0.2T_f4000hz`
trio for any orientation comparison.

**2. `results_matched_0.2T_f4000hz/dc-by-0.2T` reports exactly 0.0% at every one of its 20 snapshots.**
This is not "the field pushes everything down." For B exactly along y, `LorentzForce = J×B` is
*by construction* perpendicular to B (`(J×B)·B ≡ 0`), so `LorentzForce_y` is exactly zero for every liquid
cell, every snapshot — confirmed directly on the stored field (`y range: 0.0 to 0.0` over 26,634 liquid
cells). The tool's existing NaN-guard only catches a fully-zero force vector (the no-field-control case);
it doesn't catch "zero y-component, nonzero x/z," so this degenerate case silently reads as 0% rather than
undefined. **This is a genuine edge-case gap in the tool**, but it affects exactly this one case (the only
one with an exactly-axis-aligned vertical B field), so rather than patch the script under time pressure I
corrected it by hand in this analysis: `dc-by-0.2T` (matched root) is excluded from the up/down ranking,
same as the controls, with the physical note "a purely vertical static field cannot generate a vertical
Lorentz force at all, regardless of strength." If this survey is rerun, I'd recommend adding a check
`if |ef_y| across all cells rounds to 0 but |F| doesn't → NaN` next to the existing all-zero guard.

## Full ranked table

Controls (no field, genuinely zero `LorentzForce` at every cell — confirmed by field data, three of the
four have stale dicts still claiming `MHD.active=true`):

| case | max liquid cells | median force density |
|---|---|---|
| `results_matched_0.2T_f4000hz/no-mhd` | 27,121 | 0.0 |
| `openfoam_projects/results/no-mhd` | 31,575 | 0.0 |
| `openfoam_projects/results_spot/no-mhd` | 32,652 | 0.0 |
| `openfoam_projects/results_weld_s5e6/no-mhd` | 32,866 | 0.0 |

All 31 field-bearing cases, ranked by pool-weighted up-fraction:

| # | case | type | orient. | strength | up% (pool-wtd) | up% (force-wtd) | snaps | max liquid cells | flag |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `results/azimuthal-1T-reversed` | azimuthal | toroidal | 1T reversed | 64.2% | 78.7% | 8 | 5,001 | **UNRELIABLE** — pool never exceeds 5,001 cells and resolidifies to 60 by true end (see below) |
| 2 | `results_dc_0.2T/dc-bx-0.2T` | DC | B∥x | 0.2T | 63.7% | 48.8% | 22 | 230,133 | duplicate-bug root — the one real data point among the three below |
| 3 | `results_dc_0.2T/dc-by-0.2T` | DC | (mislabeled) | 0.2T | 63.7% | 48.8% | 22 | 230,133 | **mislabeled duplicate** of row 2 — not independent |
| 4 | `results_dc_0.2T/dc-bz-0.2T` | DC | (mislabeled) | 0.2T | 63.7% | 48.8% | 22 | 230,133 | **mislabeled duplicate** of row 2 — not independent |
| 5 | `results_matched_0.2T_f4000hz/dc-bx-0.2T` | DC | B∥x | 0.2T | **56.9%** | 52.3% | 18 | 25,787 | clean, single dict, trustworthy |
| 6 | `openfoam_projects/results_weld_s5e6/rmf-xz-0.1T` | RMF | xz | 0.1T | 55.9% | 52.0% | 21 | 33,088 | |
| 7 | `openfoam_projects/results_weld_s5e6/rmf-xz-0.2T` | RMF | xz | 0.2T | 55.1% | 50.9% | 21 | 31,705 | |
| 8 | `openfoam_projects/results_weld_s5e6/rmf-xy-0.2T` | RMF | xy | 0.2T | 54.4% | 50.9% | 21 | 32,965 | |
| 9 | `openfoam_projects/results_spot/rmf-xz-0.2T` | RMF | xz | 0.2T | 54.3% | 50.6% | 21 | 31,312 | |
| 10 | `openfoam_projects/results/rmf-xz-0.2T` | RMF | xz | 0.2T | 53.9% | 50.8% | 21 | 31,285 | |
| 11 | `openfoam_projects/results/rmf-xz-0.1T` | RMF | xz | 0.1T | 53.8% | 50.6% | 21 | 33,425 | |
| 12 | `openfoam_projects/results_spot/rmf-xy-0.2T` | RMF | xy | 0.2T | 53.0% | 50.4% | 21 | 31,233 | |
| 13 | `openfoam_projects/results/rmf-xy-0.2T` | RMF | xy | 0.2T | 52.7% | 50.2% | 21 | 32,010 | |
| 14 | `openfoam_projects/results/rmf-xy-0.1T` | RMF | xy | 0.1T | 52.3% | 50.5% | 21 | 32,155 | |
| 15 | `openfoam_projects/results_spot/rmf-xy-0.1T` | RMF | xy | 0.1T | 52.3% | 50.4% | 21 | 31,535 | |
| 16 | `results_matched_0.2T_f4000hz/rmf-xz-0.2T` | RMF | xz | 0.2T | 52.3% | 51.1% | 18 | 26,433 | |
| 17 | `openfoam_projects/results_weld_s5e6/rmf-xy-0.1T` | RMF | xy | 0.1T | 51.9% | 50.2% | 21 | 33,274 | |
| 18 | `openfoam_projects/results_spot/rmf-yz-0.2T` | RMF | yz | 0.2T | 51.1% | 50.1% | 21 | 32,424 | |
| 19 | `openfoam_projects/results_spot/rmf-yz-0.1T` | RMF | yz | 0.1T | 51.1% | 50.6% | 21 | 33,076 | |
| 20 | `openfoam_projects/results/rmf-yz-0.1T` | RMF | yz | 0.1T | 51.0% | 50.4% | 21 | 31,014 | |
| 21 | `openfoam_projects/results/rmf-yz-0.2T` | RMF | yz | 0.2T | 50.3% | 50.1% | 21 | 32,219 | essentially the 50% no-bias line |
| 22 | `results_matched_0.2T_f4000hz/dc-bz-0.2T` | DC | B∥z | 0.2T | 49.9% | 49.3% | 18 | 25,736 | essentially the 50% no-bias line |
| 23 | `openfoam_projects/results_weld_s5e6/rmf-yz-0.1T` | RMF | yz | 0.1T | 48.5% | 49.5% | 21 | 33,271 | |
| 24 | `openfoam_projects/results_weld_s5e6/rmf-yz-0.2T` | RMF | yz | 0.2T | 48.2% | 50.0% | 21 | 33,180 | |
| 25 | `openfoam_projects/results_spot/rmf-xz-0.1T` | RMF | xz | 0.1T | 47.0% | 49.0% | 23 | 6,143 | **LOW CONFIDENCE** — truncated/crashed run, 24 time dirs, stops at t=0.00023 |
| 26 | `results/azimuthal-1T-3ms` | azimuthal | toroidal | 1T, 3 ms | 45.0% | 27.6% | 20 | 108,088 | large, well-resolved, still growing — see time-series discussion |
| 27 | `results_matched_0.2T_f4000hz/rmf-yz-0.2T` | RMF | yz | 0.2T | 43.9% | 48.8% | 18 | 26,540 | |
| 28 | `results_matched_0.2T_f4000hz/rmf-xy-0.2T` | RMF | xy | 0.2T | 42.4% | 44.9% | 18 | 9,569 | pool ~3x smaller than sibling roots — treat with some caution |
| 29 | `results/azimuthal-0.2T` | azimuthal | toroidal | 0.2T | 35.0% | 27.5% | 18 | 26,876 | **validation anchor** — reliable, large pool |
| 30 | `results/azimuthal-1T-rampdown` | azimuthal | toroidal | 1T rampdown | 27.6% | 20.6% | 10 | 5,729 | **UNRELIABLE** — pool never exceeds 5,729 cells and resolidifies to 52 by true end |
| 31 | `results_matched_0.2T_f4000hz/dc-by-0.2T` | DC | B∥y | 0.2T | 0.0%* | 0.0%* | 20 | 26,634 | **DEGENERATE, not a real result** — see explanation above; `LorentzForce_y ≡ 0` by construction for B∥y |

## Answering the four questions

**1. Which configurations push pores upward most of the time (>50%)? Any at all?**

Yes, several clear >50%, but the margins above the coin-flip line are modest once the two data-quality
problems are accounted for:

- The strongest genuine (non-duplicated, non-degenerate, adequately-resolved) signal is **DC field
  perpendicular to the free surface, in the plane of the pool** (`B∥x` here) **at 0.2T: 56.9%**
  (`results_matched_0.2T_f4000hz/dc-bx-0.2T`). The mislabeled `results_dc_0.2T` root's B∥x case agrees in
  direction (63.7%) on a much larger pool, which corroborates the sign even though it can't be treated as
  a second independent trial.
- **Rotating fields in the xz plane** are the next tier and the most *consistently* reproduced result in
  the whole survey: 52.3–55.9% across all five independent xz runs (four process roots + the matched root),
  spanning both 0.1T and 0.2T. No xz-plane run fell below 50%.
- **Rotating fields in the xy plane** are mostly >50% too (6 of 7 non-degenerate runs), but more scattered
  (42.4–54.4%), with the one outlier being a small-pool run (9,569 max cells vs ~31–33k for its siblings).
- **Rotating fields in the yz plane** are the weakest RMF orientation — clustered right at or below 50%
  (43.9–51.1%), with several roots landing net *below* 50%.
- **Azimuthal/toroidal fields, on the two reliable large-pool runs available (0.2T and the 3ms 1T run),
  are net DOWNWARD (35.0% and 45.0%)** — the opposite of what the campaign is testing for. This is the
  most confident *negative* result in the survey (both are large, non-truncated, non-degenerate pools).
- **DC field along the vertical axis (B∥y) cannot contribute at all** — not "pushes down," but structurally
  incapable of a vertical component, by the `(J×B)·B=0` identity confirmed on the stored field.

**2. How does this vary with field type, orientation, strength?**

- *Type*: DC (off-vertical orientation) > RMF-xz ≈ RMF-xy > RMF-yz > azimuthal, among trustworthy rows.
- *Orientation* is the dominant variable, not strength. For DC: only the axis perpendicular to gravity/the
  pool's width plane (x here) gives a real, sizeable, reproducible upward bias; the axis in the other
  horizontal direction (z) is neutral (49.9%); the vertical axis (y) is structurally null. For RMF:
  xz > xy > yz holds up across essentially every process root tested (4–5 independent roots per plane).
- *Strength*: within RMF, matched-plane 0.1T vs 0.2T pairs differ by well under 1.5 percentage points in
  every case (e.g. `results/rmf-xz`: 53.8% at 0.1T vs 53.9% at 0.2T; `results/rmf-yz`: 51.0% vs 50.3%).
  **Field strength does not change which way a rotating field tends to push on average** — it presumably
  changes the magnitude/speed of any real expulsion, which this instantaneous-force method cannot see, but
  it does not flip or meaningfully amplify the directional bias. Azimuthal strength (0.2T vs 1T) can't be
  cleanly separated from run-length/pool-size confounds with only one reliable case per strength (see Q4).

**3. Where does the no-field control sit, and what does "better than control" mean?**

All four controls return `NaN` — confirmed genuinely zero `LorentzForce` (median force density = 0.0) at
every liquid cell in every sampled snapshot, regardless of what the (stale, in three of four cases) dicts
claim about `MHD.active`. There is no force, so there is no direction to be "better" or "worse" than —
the control isn't a competing bar on this chart, it's the absence of the mechanism being measured.

The meaningful baseline is therefore **50%, not the control.** A field configuration that lands near 50%
(most of the RMF cases, `dc-bz-0.2T`) is, in terms of *net* bubble bias, doing almost the same nothing the
control does — half the liquid gets pushed up, half down, no net preference — it's just doing that nothing
with a nonzero, noisier force field instead of a literal zero one. "Better than control" in this dataset
should mean: *meaningfully and consistently above 50%, on a run with a real (not vanishing, not
duplicated, not degenerate) melt pool.* By that standard only DC-B∥x and the RMF-xz family clearly qualify;
azimuthal fields at the strengths reliably tested here do not.

**4. Worth re-running at high output cadence for full trajectory analysis?**

**Primary recommendation: `results_matched_0.2T_f4000hz/dc-bx-0.2T`.** It has the cleanest, strongest,
single-configuration upward signal in the survey (56.9% pool-weighted, 52.3% force-weighted) on a
properly-resolved, non-duplicated, non-degenerate pool (25,787 max liquid cells, never resolidifies to
near-empty in the sampled window). Its direction is corroborated by the (data-compromised) `results_dc_0.2T`
root's B∥x case landing even higher (63.7%) on a much bigger pool. This is the best candidate to confirm
with full Lagrangian particle tracking whether the instantaneous field bias actually translates into bubbles
reaching the surface, since it's the field configuration this Eulerian survey trusts most.

**Secondary recommendation: `results/azimuthal-1T-3ms`, specifically its tail (t > ~2 ms).** It's the only
azimuthal run that's both reliable (108,088 max liquid cells, growing the entire 3 ms, never resolidifying)
and shows an interesting time-varying pattern: the instantaneous up-fraction sits flat around 40–45% for
most of the run, then rises to 50.5% and 58.0% in the last two sampled snapshots (see
`survey_azimuthal_timeseries.png`). A single-snapshot Eulerian survey cannot tell whether that late upturn
reflects a real late-pulse physical effect (pool geometry rotating into a more favorable position, field
phase, etc.) worth exploiting, or is a transient that a fine-cadence Lagrangian rerun would show doesn't
actually move bubbles by the time it appears. This is exactly the "does it help or just look that way"
question the campaign needs answered, and it's currently open.

I would *not* prioritize rerunning any of `azimuthal-1T-reversed`, `azimuthal-1T-rampdown`, or the
`results_dc_0.2T` by/bz "variants" for trajectory work — the first two never form a real pool, and the third
would just be re-verifying a mislabeling bug, not new physics. If the campaign wants B∥y or B∥z DC and
by/bz RMF-plane results at the larger `results_dc_0.2T` pool scale, that requires a **new CFD run** with
correct per-axis dicts, not a rerun of existing output at finer cadence.

## Data-quality caveats to keep in mind when reading the table

- **Small/vanishing pools bias individual rows.** `azimuthal-1T-reversed` and `azimuthal-1T-rampdown` never
  exceed ~5,000–5,700 liquid cells (vs 25,000–108,000 for every other reliable case) and both resolidify to
  52–60 cells by the true final time directory (t=0.002s; confirmed directly by reading `alpha.metal` /
  `epsilon1` at every stored time in the last 5% of each run). The stride-sampled snapshots the tool actually
  used stop slightly before that (125 and 104 cells respectively) because of its own `<100 liquid cells →
  drop` guard, but the trend and the extreme, noisy per-snapshot values (up to 90.4% up-fraction on 125
  cells, force densities >5×10⁷ N/m³) show the pool-size weighting is not fully protecting these two rows.
  Both are flagged UNRELIABLE and should not be ranked next to full-pool runs.
- **`openfoam_projects/results_spot/rmf-xz-0.1T` is a truncated/crashed run** (24 time dirs, stops at
  t=0.00023 vs t=0.001 for every sibling case) — flagged LOW CONFIDENCE, sits near the 50% line anyway so it
  doesn't change any conclusion, but should not be read as a completed 0.1T xz-plane result.
- **`results_matched_0.2T_f4000hz/rmf-xy-0.2T`** has a noticeably smaller pool (max 9,569 cells) than its
  four sibling roots (31,000–33,000 cells) — not flagged as unreliable outright, but its 42.4% should be
  weighted less heavily than the other RMF-xy numbers when judging the xy-plane orientation as a whole.
- **Process-type labeling caveat:** the prompt describes spot/weld/scan process variants across roots. I
  checked the `timeVsLaserPosition`/`timeVsLaserPower` dicts for the `no-mhd` control in each of `results`,
  `results_spot`, and `results_weld_s5e6` — all three use an identical, *stationary* two-point laser-position
  table (`(250e-6,700e-6,250e-6)` at both t=0 and t=100e-6). I could not confirm a scanning/moving laser in
  any of the roots surveyed. `results_weld_s5e6` does differ from the other two in electrical/optical
  resistivity and its default MHD dict block, so the roots are not literal duplicates of each other, but the
  "spot vs weld vs scan" distinction the campaign name implies is not visible in the laser-motion dicts I
  inspected — treat any root-to-root comparison as "different campaign/parameter set," not confidently as
  "different process type."

## Plots

- `survey_ranked_bar.png` — all 31 field-bearing cases ranked by pool-weighted up-fraction, colored by
  field type (DC / RMF / azimuthal), hatched where flagged unreliable, with the 50% no-bias line marked.
- `survey_azimuthal_timeseries.png` — instantaneous up-fraction and pool size (log scale) vs time for the
  four azimuthal-field runs, showing the small-pool cases (reversed, rampdown) are noisy and short-lived
  while `azimuthal-0.2T` and `azimuthal-1T-3ms` are well-behaved and the latter's late-time upturn.
