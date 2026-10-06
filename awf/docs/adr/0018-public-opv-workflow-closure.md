# ADR 0018: Publish the complete OPV workflow closure

Status: accepted for `gromacs-organic-blend-zc@1.0.0-rc.6` (2026-10-05).

## Context

The sanitized OPV Artifact in Workflow rc.4 removed private paths and Gaussian
files from the data package. A public-package audit then found that the locked
Tool and production Recipe still contained internal filesystem provenance. A
public Workflow must be assessed as its complete locked dependency closure,
not as the data Artifact alone.

Workflow rc.5 selected neutral successors for those two Components, but its
engineering-gate text retained the internal Artifact's old checksum count.
Because registered Releases are immutable, that wording is corrected in a new
Workflow Release rather than edited in place.

## Decision

`gromacs-organic-blend-zc@1.0.0-rc.6` locks the public Artifact rc.2, Tool rc.2
and production Recipe rc.2. The remaining locked Recipe and Skills were scanned
as part of the same closure. Public packaging excludes generated caches and
fails if it finds private paths, private contributor labels, Gaussian input or
output files, or Python bytecode.

The public engineering gate records 218 Artifact file checksums. Release rc.5
is retained exactly as registered; rc.6 becomes the Registry current pointer.

## Consequences

The GitHub package can be reproduced from one exact lock and verified without
access to the original project tree. The sanitized package preserves 91
fragment records and 10 molecule packages. Scientific-data rights remain
separate from the MIT license for code and metadata, as stated in the bundled
data notice.
