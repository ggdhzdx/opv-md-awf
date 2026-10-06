# GROMACS organic-blend production recipe

This recipe defines a matched, replicated workflow for atomistic organic-blend
morphology simulations.  Molecular identities, force-field files, coordinates,
machine paths, and trajectories are project data and are not bundled here.

## Required project adapter

The project must provide a machine-readable manifest containing:

- two or more named blend systems;
- donor and acceptor coordinate/topology paths plus SHA-256 hashes;
- molecule counts, coordinate order, topology order, initial-box size, and
  Packmol tolerance;
- at least three independent replica seeds, with the same Packmol seed reused
  across matched chemistries within a replica and no coordinate/checkpoint
  ancestry shared between replicas;
- explicit GROMACS build, accelerator, host, and resource layout;
- the stage protocol and quantitative acceptance gates.

## Scientific and technical gates

1. Confirm force-field provenance, net charge, atom count/order, connectivity,
   and coordinate/topology identity before packing.
2. Treat static `grompp`, energy-minimization startup, and runtime startup checks
   as mandatory technical validation. They are not a smaller scientific pilot.
3. Generate each replica independently. A matched seed controls the comparison;
   it does not make different chemistries coordinate-isomorphic.
4. Run energy minimization, staged compression, 300 K/1 bar relaxation, three
   anneal cycles, and a 300 K/1 bar hold. The default schedule is in
   `protocol.json`; a project may override it only with explicit scientific
   authorization and a new manifest hash.
5. Advance a cohort together through density/QC gates. Never silently drop a
   failed chemistry or replica and continue the comparison.
6. Resume only from a checkpoint whose hashes and ancestry match the immutable
   run manifest. Partial checkpoints, mismatched TPRs, sparse text summaries,
   and coordinate/topology drift fail closed.
7. Report results as fixed-protocol, quenched-morphology observations. Do not
   claim equilibrium morphology or causal device performance from this workflow
   alone.

## Acceptance evidence

Retain the frozen input manifest, environment receipt, Packmol inputs/logs,
per-stage `grompp` and `mdrun` logs, checkpoint ancestry, energy files, density
and stability metrics, failure ledger, and final analysis input hashes.

