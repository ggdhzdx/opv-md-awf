---
name: run-gromacs-organic-blend
description: Prepare, run, resume, and audit matched replicated GROMACS morphology simulations for atomistic organic donor-acceptor blends. Use when a project needs Packmol initialization, staged compression/annealing/hold production, checkpoint-safe continuation, or comparison across related blend chemistries; do not use for force-field-free coarse-grained models or single-crystal lattice dynamics.
---

# Run GROMACS Organic Blends

Build a project-specific adapter around the locked AWF recipe. Keep all
molecular inputs and scientific results in the project.

## Before execution

1. Read the project authorization and the applicable cluster/program skills.
2. Read the locked recipe `RECIPE.md` and `protocol.json` completely.
3. Inspect existing validated donor/acceptor topology resources before creating
   new parameters. Reuse only graph- and provenance-compatible fragments.
4. Require explicit evidence for atom order, bond graph, net charge, molecule
   counts, topology order, and force-field method. File integrity alone is not
   chemical compatibility.
5. Query the program run-guide index for the exact GROMACS build, hardware,
   backend, task class, and system scale. If no compatible guide exists, record
   a baseline reset; do not declare an unbenchmarked layout optimal.

## Execution contract

- Freeze the project manifest and all input hashes before Packmol.
- Use at least three independent replicas. Reuse a Packmol seed only across
  matched chemistries within one replica; never share coordinates or
  checkpoints across replicas.
- Preserve explicit molecule and topology order and audit every packed molecule.
- Perform mandatory `grompp`, minimization-start, and production-start checks.
  These are technical correctness gates, not a reduced scientific pilot.
- Advance all matched systems through the same stage schedule and cohort gate.
- Keep stage markers, TPR/MDP/input hashes, checkpoint ancestry, exit status,
  energy files, and logs. Resume only on an exact ancestry match.
- On a new error, leave running jobs unchanged, collect evidence, and report
  before changing the method, input, parameters, or resource layout.

## Analysis and claims

Bind analysis to hashed TPR/trajectory pairs and report replica-level values
before aggregate comparisons. Preserve uncertainty and paired design. Describe
the output as conditioned quenched morphology unless independent evidence
supports equilibrium. Keep correlations with device performance separate from
causal claims.

## Candidate boundary

This rc.1 release captures workflow decisions and gates. It is not adoptable
until a second project adapter, negative regressions, and a real replicated
production run verify the locked behavior.

