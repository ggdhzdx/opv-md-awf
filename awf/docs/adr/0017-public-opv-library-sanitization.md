# ADR 0017: Publish a sanitized OPV library successor

Status: accepted for `gromacs-organic-blend-zc@1.0.0-rc.4` (2026-10-05).

## Context

The internal OPV library Artifact made `opv-md` portable, but it retained
Gaussian input/output provenance, private absolute paths and source-set labels
derived from personal project directories. Those bytes were useful internally
but inappropriate for a public GitHub distribution.

## Decision

`artifact.opv-fragment-library@1.0.0-rc.2` is a new immutable public successor.
It removes the Gaussian directory, replaces private paths with package-relative
or redacted provenance references, maps person-identifying source labels to
neutral source ids, and regenerates every catalog and checksum. It retains all
91 indexed fragment records, 10 molecule packages and the coordinate/topology
files required by the library Tool.

The AWF code license does not automatically grant rights in the scientific data;
the Artifact carries a separate `DATA-NOTICE.md`. A later explicit owner decision
may publish a separately licensed data successor without changing this Release.

## Consequences

`gromacs-organic-blend-zc@1.0.0-rc.4` can be distributed without Gaussian log
license text or private filesystem/person labels. The internal rc.1 Artifact and
Workflow rc.3 remain reproducible. Sanitization preserves reference integrity,
not chemical validity or a general redistribution-rights determination.
