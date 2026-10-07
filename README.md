# opv-md AWF public distribution

Portable public distribution of `gromacs-organic-blend-zc@1.0.0-rc.7`, alias `opv-md`. It prepares, audits, resumes, and analyzes replicated atomistic organic-blend morphology simulations.

This successor fixes the protocol-only analysis closure in rc.6. It now contains runnable BulkGraph and phase-separation programs. The historical `candidate` label remains audit metadata and does not gate task use.

## Included

- exact Workflow descriptor and version lock;
- two Skills, two Recipes, the OPV library Tool, and the fragment-library Artifact;
- 91 indexed fragment records, 10 molecule packages, and 218 internal checksums;
- `opv-bulkgraph`: NumPy NPZ/JSON contact and stacking analysis;
- `opv-phase-separation`: optional MDAnalysis-backed GROMACS TPR/XTC CIC/FFT analysis;
- synthetic numerical tests, scientific boundaries, and standalone package verification.

## Install and run analysis

```bash
cd awf/components/recipes/opv-morphology-analysis/releases/1.0.0-rc.2/source
python -m pip install .
opv-bulkgraph --manifest bulkgraph-input.json --output bulkgraph-result.json

# Direct GROMACS trajectory support
python -m pip install '.[gromacs]'
opv-phase-separation --mode full_run --tpr production.tpr --xtc production.xtc   --out-dir phase-results --donor-selection 'resname DON' --acceptor-selection 'resname ACC'
```

See the locked Recipe for the BulkGraph manifest contract and analysis gates.

## Verify

```bash
python verify.py
```

The verifier checks every distributed file, the exact lock, both executable entry points, Python/JSON syntax, public-content exclusions, and the complete fragment-library inventory.

Private filesystem provenance, fixed project molecule identities, target trajectories, checkpoints, results, credentials, and external programs are excluded. GROMACS, Packmol, ztop, and project inputs remain external.

File identity and synthetic numerical regression do not establish force-field quality, equilibrium morphology, electronic coupling, or device performance.
