---
name: analyze-opv-morphology
description: Execute and audit atomistic organic photovoltaic (OPV) donor-acceptor morphology analysis, integrating BulkGraph contact/stacking networks and nanoscale phase-separation / composition heterogeneity structure factors. Use when a project requires molecular contact graphs, SVD planar stacking, Bondi vdW cofacial overlap area, or periodic Cloud-in-Cell structure factors and coarse-grained composition variance; do not use for coarse-grained implicit solvent models or single-crystal lattice dynamics.
---

# Analyze OPV Morphology

Build a project-specific adapter around the locked AWF recipe. Keep all
molecular inputs, trajectories, and scientific results in the project.

## Before execution

1. Read the project authorization, method confirmations, and applicable cluster run guides.
2. Read the locked recipe `RECIPE.md` and `protocol.json` completely.
3. Verify that the simulation trajectory originates from an audited production schedule (e.g. `gromacs-organic-blend-zc@1.0.0-rc.1`).
4. Require two-phase read-only SHA-256 digest verification on source `.tpr` and `.xtc` files against `_COMPLETE.json` before and after analysis.
5. Query the program run-guide index (`~/.codex/program-run-guides/index.yaml`) for the matching execution profile on the target hardware (e.g. `node01_zen3` with `OPENBLAS_NUM_THREADS=1`). Do not launch unbenchmarked threading layouts.

## Execution contract

- Parallel Analysis Modules:
  * Module A (BulkGraph):
    - Compute molecular contact graph $G_{\text{contact}}$ on heavy atoms of conjugated fragments ($\le 0.45\text{ nm}$).
    - Perform SVD planar fitting on conjugated cores (`ADN_core`, `AAC_left`, `AAC_right`).
    - Classify primary stacking: $\theta \le 30.0^\circ$, normal projection $\text{proj} \ge 0.25\text{ nm}$ with antiparallel sign alignment, $d_{\min} \le 0.36\text{ nm}$.
    - Classify sensitivity sets ($d_{\min} \le 0.40\text{ nm}$, $\theta \le 25.0^\circ$, $\theta \le 35.0^\circ$).
    - Evaluate Bondi vdW projected cofacial overlap area strictly on primary stacked pairs with adaptive grid convergence ($h \to h/2 \to h/4$).
  * Module B (Phase Separation & Composition Heterogeneity):
    - Assign mass to periodic 3D grids using trilinear Cloud-in-Cell (CIC) interpolation.
    - Evaluate local mass fraction $\phi_{\text{NFA}}(\mathbf{r})$ and spatial fluctuation $\delta\phi(\mathbf{r})$ centered by the spatial grid mean.
    - Compute concentration structure factor $S_{cc}(q)$ via 3D FFT and spherical shell averaging aligned on contiguous common integer shells.
    - Enforce the characteristic scale resolution rules: require an interior peak ($b^* > b_{\text{qmin}}$), discrete curvature, and multi-grid robustness; classify monotonically decaying lowest-$q$ maxima as `box_scale_limited`, and otherwise report `no_resolved_characteristic_scale` ($q^* = \text{null}$, $d = \text{null}$).
    - Compute Gaussian-filtered coarse-grained composition variance $\sigma_\phi^2(\sigma)$ at $\sigma \in \{1.0, 2.0, 3.0\}\text{ nm}$.

## Analysis and claims

- Report replica-level time-averaged values before cohort contrasts.
- Preserve paired replica design ($\Delta_{\text{CF3}-\text{CH3}}$, $\Delta_{\text{Ph}-\text{CH3}}$, $\Delta_{\text{Ph}-\text{CF3}}$) and document linear dependence.
- Characterize results strictly as fixed-protocol conditioned quenched morphology and nanoscale composition heterogeneity.
- Barred inferences:
  * NO claim of thermodynamic equilibrium or Flory-Huggins $\chi$ parameters.
  * NO claim of macroscopic phase separation (finite ~10 nm box over 10 ns cannot demonstrate macroscopic phase separation).
  * NO causal claims linking contact/stacking/heterogeneity metrics to device solar cell performance ($V_{\text{oc}}, J_{\text{sc}}, \text{FF}, \text{PCE}$).

## Candidate boundary

This rc.1 release captures workflow decisions and mathematical gates. It is not adoptable
until multi-system replicated production acceptance verifies the locked behavior.
