# ADR 0016: Bundle small reusable scientific Artifacts

Status: accepted for `gromacs-organic-blend-zc@1.0.0-rc.3` (2026-10-05).

## Context

`opv-md` selected a library Tool but kept the 29.8 MiB OPV fragment and
molecule library outside AWF. A transported Workflow therefore contained the
search and validation code without the 91 fragment records and 10 molecule
packages it was meant to inspect. The Workflow was reproducible only on the
original project filesystem.

## Decision

AWF may bundle a small, redistributable, source-preserving scientific resource
as an immutable Artifact source snapshot when that resource is required for a
portable Workflow. The Artifact must record provenance, byte count and tree
digest, retain its internal checksum manifest, and state its scientific limits.
Large runtimes, model payloads, project inputs, trajectories, checkpoints and
results remain external.

`artifact.opv-fragment-library@1.0.0-rc.1` freezes the verified 2026-09-15
library. `gromacs-organic-blend-zc@1.0.0-rc.3` selects it in `versions.lock`.

## Consequences

The current `opv-md` dependency closure is self-contained for fragment lookup,
reference checks and reuse planning. Core bundles grow by about 30 MiB. Library
integrity does not establish attachment-site correctness, charge, force-field
compatibility, single-molecule stability, morphology or device performance.
The original project library and all earlier AWF Releases remain unchanged.
