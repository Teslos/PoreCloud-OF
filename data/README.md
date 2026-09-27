# Curated inputs

Small, hand-checked inputs the analysis depends on. Not simulation output —
`.gitignore` carries an explicit negation so the blanket `*.csv` rule does not
swallow them.

| directory | what it is |
|---|---|
| `calibration/` | per-case injection rates. `rate` is a nucleation rate per unit **liquid** volume per second, so the resulting count is `rate · V_melt · T` and it must be calibrated per case. An uncalibrated value once produced bubbles occupying 473% of the melt volume, which the solver does not object to because the coupling is one-way. |
| `survey/` | the 35-case Eulerian force-direction survey (`poreCloudFieldSurvey.py`) |
| `poretracker-summaries/` | resolved-VoF pore summaries from `results_matched_0.2T_f4000hz`. The source of **both** the injection calibration and the radial seeding weights. |
| `invalid-runs/` | notes on runs kept as a record but whose numbers must not be quoted |
