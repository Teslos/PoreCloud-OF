# Production case setups

The ten `poreCloud-*` directories are the configurations behind the expulsion
campaign, carried as **dictionaries only**. Each holds `system/` and `constant/`
and deliberately no mesh, no time directories and no `processor*` — the runs they
came from total ~440 GB, and none of it is needed to reproduce them.

`initial-common/` is the cold-start field set. All seven matched-family cases
used byte-identical initial fields, and every `internalField` in them is
`uniform`, so the one copy is mesh-independent and serves every case here.

## Running one

    cd <case>
    blockMesh                      # regenerates constant/polyMesh
    cp -r ../initial-common 0
    setFields
    decomposePar && mpirun -np 18 laserbeamFoam -parallel

## One thing to know about the numbers in docs/

The published runs were **restarts from t = 1 ms** of a parent run, so the melt
pool was already developed when the cloud started injecting; their
`controlDict` says `startFrom latestTime` with `endTime 0.002`. Starting cold
from `initial-common` instead runs 0 → 2 ms and spends the first millisecond
developing the pool, which roughly doubles the wall time and does not reproduce
the parent's exact turbulent state. For a like-for-like comparison against
`docs/specs.md`, restart from a parent at t = 1 ms.

## The suffixes

| suffix | meaning |
|---|---|
| `-cal` | injection rate calibrated for that case (see `../data/calibration/`) |
| `-seeded` | births drawn from the measured radial distribution, not uniformly |
| `-fix` | re-run after a correctness fix; supersedes the unsuffixed run |
| `_INVALID-*` | kept as a record, **not** to be quoted — reason in the name |

`poreCloud-ac-by-400Hz` is superseded by `poreCloud-ac-by-400Hz-fix`; see
`../data/invalid-runs/` and `docs/specs.md` §5.13 for why.
