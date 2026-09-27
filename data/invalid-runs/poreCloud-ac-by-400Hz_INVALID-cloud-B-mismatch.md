# INVALID — cloud and solver used different B fields

Completed t = 1 -> 2 ms, 376 parcels, exit 0. The **melt flow is correct**; every
**cloud force and trajectory is not**.

`poreCloudMagneticField.C` had no `alternating` branch. An AC case sets both
`frequency` and `rotationPlane`, so `buildUniform()` fell through to the RMF
branch and gave the cloud

    B = 0.2 (cos wt, sin wt, 0)      # rotating in the xy plane

while `UEqn.H` applied to the momentum equation

    B = 0.2 (0, cos wt, 0)           # oscillating along y

The tell: this run reports max|F_y| = 3.0e-6 N, the same order as
max|F_x| = 4.7e-6 N. With B along y the vertical exclusion force is
*identically* zero, because (J x B)_y = J_z B_x - J_x B_z and both B_x and B_z
vanish. A nonzero F_y was therefore impossible for the field that was actually
applied, and the run appeared to produce vertical force where there was none.

Kept as the record of the discrepancy. Do not quote its numbers.
Fixed in PoreCloud-OF by adding the `alternating` branch to `buildUniform()`.
