# Provenance

This public successor turns the rc.1 mathematical contract into executable,
project-neutral Python modules. The algorithms were generalized from audited
atomistic donor/acceptor morphology work, with project paths, molecule names,
trajectories, results, and fixed chemistry schemas excluded.

The executable closure contains a NumPy BulkGraph implementation, periodic CIC
mass assignment, concentration structure factors, multi-grid peak gates,
Gaussian-filtered composition variance, a portable NPZ/JSON BulkGraph adapter,
and an optional MDAnalysis adapter for GROMACS TPR/XTC trajectories.

No scientific input data or third-party runtime is bundled.
