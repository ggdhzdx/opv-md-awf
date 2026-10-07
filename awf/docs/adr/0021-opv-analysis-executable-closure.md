# ADR 0021: Make OPV morphology analysis executable

Status: accepted for `gromacs-organic-blend-zc@1.0.0-rc.7` (2026-10-07).

## Context

The public rc.6 closure contained the OPV morphology Skill, Recipe, protocol,
and mathematical gates, but it contained no executable analysis program. A
recipient could inspect the intended BulkGraph and phase-separation methods but
could not run them from the repository.

## Decision

Recipe `opv-morphology-analysis@1.0.0-rc.2` adds two public entry points:

- `opv-bulkgraph`, a NumPy NPZ/JSON adapter for contact and stacking graphs,
  SVD plane fitting, sensitivity gates, and optional projected overlap area;
- `opv-phase-separation`, an optional MDAnalysis-backed GROMACS TPR/XTC adapter
  for periodic CIC composition fields, spherical FFT structure factors,
  characteristic-scale gates, and Gaussian-filtered variance.

The public implementation removes project paths, fixed molecule identities,
trajectories, results, and private chemistry schemas. BulkGraph currently fails
closed for non-orthorhombic cells. Skill rc.2 names the executable entry points,
and Workflow rc.7 locks both successor components. Earlier Releases remain
immutable.

## Consequences

The public repository now provides runnable analysis rather than protocol-only
guidance. Synthetic numerical regression establishes code-path behavior and
failure gates. It does not validate a project's atom selections, sampling,
force field, ensemble convergence, or scientific interpretation.
