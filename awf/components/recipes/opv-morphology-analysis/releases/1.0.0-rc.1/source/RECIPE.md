# OPV morphology analysis recipe

This recipe defines a matched, replicated workflow for atomistic organic photovoltaic
(OPV) donor-acceptor blend morphology analysis. Molecular identities, trajectories,
force-field files, coordinates, and machine paths are project data and remain in the project.

## Required project adapter

The project must provide a machine-readable manifest containing:

- two or more named blend systems (e.g. `D18_A-CH3`, `D18_A-CF3`, `D18_A-Ph`);
- at least three independent replicas (`r01`, `r02`, `r03`);
- trajectory paths with two-phase SHA-256 digest validation against production manifests;
- analysis window definition: 5.0–10.0 ns (5000.0–10000.0 ps) with 20.0 ps cadence (251 frames);
- atom selection definitions for conjugated cores (`ADN_core`, `AAC_left`, `AAC_right`) and component species;
- execution environment profile matching the calibrated program run guide (`OPENBLAS_NUM_THREADS=1`).

## Parallel analysis modules

1. **Module A: BulkGraph Contact & Stacking Network**:
   - Molecular Contact Graph ($G_{\text{contact}}$): Intermolecular heavy-atom contact cutoff $\le 0.45\text{ nm}$ between conjugated fragments.
   - SVD Plane Fitting: Normal vectors fit to heavy atoms of conjugated rings; antiparallel sign alignment ($\mathbf{n}_b' = \operatorname{sgn}(\mathbf{n}_a \cdot \mathbf{n}_b)\mathbf{n}_b$).
   - Primary Stacking Criteria: Tilt angle $\theta \le 30.0^\circ$, normal projection $\text{proj} \ge 0.25\text{ nm}$, heavy-atom minimum distance $d_{\min} \le 0.36\text{ nm}$.
   - Sensitivity Sets: Distance $d_{\min} \le 0.40\text{ nm}$; angle sensitivities at $25.0^\circ$ and $35.0^\circ$ maintaining subset inclusions ($\text{angle\_025} \subseteq \text{primary} \subseteq \text{angle\_035}$).
   - Projected Cofacial Overlap Area: Evaluated strictly on primary stacked pairs using Bondi vdW disk unions on the mutual mean plane with cell-centered adaptive grid refinement ($h \to h/2 \to h/4$).
   - Fail-closed condition: If area fails to converge, emit `status = "BLOCKED_NUMERICAL_CONVERGENCE"` with `area_nm2 = None`.

2. **Module B: Phase Separation & Nanoscale Composition Heterogeneity**:
   - Density Fields: Periodic Cloud-in-Cell (CIC) trilinear mass assignment on regular grids ($\Delta x \approx 0.20\text{ nm}$ primary, $0.25\text{ nm}$ secondary).
   - Local Mass Fraction: $\phi_{\text{NFA}}(\mathbf{r}) = \rho_{\text{NFA}}(\mathbf{r}) / (\rho_{\text{NFA}}(\mathbf{r}) + \rho_{\text{D18}}(\mathbf{r}))$.
   - Spatial Centering: Fluctuation field $\delta\phi(\mathbf{r}) = \phi_{\text{NFA}}(\mathbf{r}) - \langle \phi_{\text{NFA}} \rangle_{\text{spatial}}$ guarantees zero DC mode.
   - Structure Factor $S_{cc}(q)$: 3D FFT power spectral density spherically averaged into contiguous integer shells aligned across frames from the lowest common shell $b_{\text{qmin}}$.
   - Characteristic Scale Rules:
     * Interior peak required ($b^* > b_{\text{qmin}}$) with discrete curvature and multi-grid consistency.
     * If lowest-$q$ shell is global maximum with monotonic decay, classify as `"box_scale_limited"` ($q^* = \text{null}$, $d = \text{null}$).
     * Otherwise classify as `"no_resolved_characteristic_scale"` ($q^* = \text{null}$, $d = \text{null}$).
   - Coarse-Grained Composition Variance: $\sigma_\phi^2(\sigma)$ at $\sigma \in \{1.0, 2.0, 3.0\}\text{ nm}$ via reciprocal Gaussian filtering.

## Scientific and technical gates

1. Single-frame evaluations (`trajectory[-1]`) are strictly prohibited; 251-frame time averages are mandatory.
2. Incomplete metric extractions must fail closed (KeyError/ValueError); partial cohort summaries are prohibited.
3. No claims of thermodynamic equilibrium, Flory-Huggins $\chi$ parameters, macroscopic phase separation, or causal device efficiency.
