# opv-md AWF public distribution

This repository is the portable public distribution of `gromacs-organic-blend-zc@1.0.0-rc.6`, alias `opv-md`. It prepares, audits, resumes and analyzes replicated atomistic organic-blend morphology simulations.

The exact Release is recorded because AWF uses a rolling current pointer. The historical `candidate` label is retained for compatibility and audit history; it does not block task-scoped use.

## Included

- exact Workflow descriptor and version lock;
- two AWF Skills, two Recipes and the OPV library Tool;
- `artifact.opv-fragment-library@1.0.0-rc.2`;
- 91 indexed fragment records, 10 molecule packages and 218 internally checksummed public files;
- neutral source identifiers, collision records and scientific-use boundaries.

Private filesystem provenance, personal source labels, quantum-chemistry input/output, target-project trajectories, checkpoints, results, credentials and external programs are excluded. GROMACS, Packmol, ztop and a project adapter must be supplied separately.

## Verify

```bash
python verify.py
```

This verifies the repository manifest, the exact Workflow lock, the public-content exclusions, and all fragment-library references and checksums.

## Use with AWF

`awf/registry.fragment.json` is merge input, not a replacement registry. Publish the listed immutable Releases through the recipient AWF installation's normal process, merge the fragment, and run that installation's verifier before use.

File and reference integrity does not establish chemical compatibility, force-field quality, equilibrium morphology, electronic coupling or device performance.
