# Provenance of OPV morphology analysis recipe

This recipe was generalized from the audited morphology analysis workflows developed for the D18/three-NFA molecular dynamics campaign (`p1-n3-d18-three-nfa-20260920a`), specifically:

1. **BulkGraph Analysis Module**:
   - Developed in `d18_nfa_md/analysis_consumer_v1/`.
   - Confirmed under Rev 5 Stacking Angle Sensitivity and Distance Sensitivity specifications.
   - Provides molecular contact graph, SVD planar fitting, primary stacking, 3-level angle sensitivity hierarchy ($\text{angle\_025} \subseteq \text{primary} \subseteq \text{angle\_035}$), and fail-closed grid overlap area.

2. **Nanoscale Phase-Separation & Composition Heterogeneity Module**:
   - Developed in `d18_nfa_md/phase_separation_analysis_v1/`.
   - Confirmed under Rev 1 Phase-Separation Protocol.
   - Provides mass-weighted periodic Cloud-in-Cell density fields, spatial fluctuation centering, 3D FFT concentration structure factor $S_{cc}(q)$ with shell alignment, explicit box-scale-limited / unresolved domain scale gates, and coarse-grained composition variance $\sigma_\phi^2(\sigma)$.

Project molecular structures, topologies, force fields, simulation trajectories, and production results are excluded from this AWF recipe.
