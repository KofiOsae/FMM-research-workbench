# Research readiness and ten-item roadmap

Status date: 2026-09-28. A feature is marked **complete** only when its solver, interface, tutorial, and automated reference checks cover the stated scope.

| # | Required improvement | Current status | Exact scope and remaining gate |
|---|---|---|---|
| 1 | Patterned-stack complex S-matrix poles | Missing | The available pole tool is restricted to uniform multilayers with frozen optical indices. A patterned outgoing-wave complex-frequency solver is still required. |
| 2 | Quasinormal-mode fields and normalization | Missing | No normalized patterned QNM field is returned. The field workspaces are driven-frequency fields or closed/bound eigenfields. |
| 3 | Radiative, absorptive, and total Q for patterned structures | Partial | Far-field linewidth fitting and uniform-stack loss-removal estimates exist. Patterned pole-based Q decomposition is not implemented. |
| 4 | Ideal-BIC identification and symmetry classification | Partial | Signed-angle, symmetry, convergence, Q-versus-angle/asymmetry, and field diagnostics are available. They do not rigorously certify an ideal BIC. |
| 5 | Polarization-vortex winding and topological charge | Partial | Closed-loop winding with complex sideband background subtraction is available as a stability diagnostic. Certification still requires radiation coefficients from an isolated patterned pole or eigenmode. |
| 6 | Full-vector 2D waveguide modes | Partial | A rectangular, isotropic dielectric cross section returns all six complex components using a vector finite-difference eigenproblem. Arbitrary geometry, tensor media, open/PML boundaries, are still required. A three-run mesh/domain comparison and numerical history are now available; mode identity and acceptance tolerance still require review. |
| 7 | Purcell factor, LDOS, beta factor, and collection efficiency | Partial | Planar multilayer electric-dipole LDOS and upper-objective collection are implemented. Patterned Green functions, normalized-mode projection, named-mode beta, and patterned collection remain missing. |
| 8 | Correlated and Bayesian uncertainty | Partial | Correlated fabrication tolerance and bounded multi-chain Metropolis inference are implemented. The likelihood is independent Gaussian and diagnostics are basic; correlated likelihoods, hierarchical models, rank-normalized diagnostics, and posterior predictive checks remain. |
| 9 | General constrained geometry optimization | Partial | One-to-five bounded variables, linear constraints, multiobjective targets, differential evolution, local polishing, and checkpointed sweeps are available. Arbitrary topology/shape optimization, gradients or adjoints, discrete fabrication rules, and robust uncertainty-aware objectives remain. |
| 10 | Thermal, nonlinear, gain, anisotropic, and magnetic physics | Partial | A homogeneous local constitutive calculator handles rotated uniaxial permittivity, thermo-optic and Kerr shifts, scalar permeability, and small-signal gain. It is not coupled self-consistently to FMM, heat flow, nonlinear iteration, gain saturation, or tensor interfaces. |

## Cross-cutting resonance-fit requirement

The original adaptive fitter tracks one selected maximum or minimum. A separate multi-feature workflow now implements the following far-field operations:

1. detect maxima and minima on a global scan with a user-selected maximum count;
2. refine every retained candidate independently and jointly fit all components with a shared constant, linear, or quadratic background;
3. compare resonance counts and Lorentzian/Fano combinations using BIC and residual RMSE;
4. report each center, linewidth, loaded Q, covariance-derived uncertainty, and the full parameter-correlation matrix;
5. flag linewidth uncertainty, spectral overlap, strong cross-component correlation, and competing models with small ΔBIC;
6. assign branches using linewidth-normalized wavelength continuity plus an optional normalized coarse field signature; and
7. label every result as a phenomenological far-field fit, explicitly excluding poles, QNMs, and radiative/absorptive Q decomposition.

Remaining work is to connect the multi-branch assignment directly to every general sweep/optimizer workflow and to add pole-residue field overlap after the patterned outgoing-mode subsystem exists.

## Interface and documentation audit

| Requirement | Status | Evidence and action required |
|---|---|---|
| Tutorial covers every available workspace | Partial | The tutorial covers the main solver distinctions and new advanced cards. It still needs a searchable parameter glossary, complete worked examples for every advanced card, and synchronized acceptance tests for all instructions. Add an illustrated band-diagram guide showing gaps, flat and steep bands, crossings and avoided crossings, mode localization, wavelength scaling, and finite-slab limits. |
| Every 2D quantitative map has axes, units, title, numeric colorbar, hover values, and working palette choice | Partial | The main scattering and several field maps meet most of this contract. Vector modes and several secondary maps do not yet share the complete renderer; some palette selectors recolor only the legend rather than the raster. |
| Professional task organization | Partial | Workspaces, calculation selectors, project save/open, source-aware materials, methods export, reference cases, and warnings are present. Long cards, inconsistent result controls, lack of searchable help, and lack of a background-job manager still impede routine use. Redesign crowded multi-series legends and the measured/simulated overlay annotation together in a later visual-consistency pass. |
| Active-app feature parity | Partial | A direct **Track several resonance branches** shortcut opens the existing sweep tracker; persistent run controls expose workspace calculations. A dedicated published polariton tracking exercise remains to be added. |
| Every parameter explains meaning, units, valid range, numerical impact, and experimental choice | Partial | Major FMM inputs are explained. New mode, Bayesian, pole, constitutive, and optimization parameters need contextual definitions, recommended ranges, and consequence warnings beside each control. |
| Literature and method provenance | Partial | Core FMM, materials, bands, GMR, energy, BIC, vector-mode, planar Green-function, Bayesian-diagnostic, and QNM sources are listed. Each solver result still needs a direct method citation in its export manifest and additional primary references for coupled thermal, nonlinear, gain, and magnetic implementations before those solvers are built. |

## Release rule

Do not describe the application as generally research-grade or as a Lumerical replacement. A result becomes suitable for a research claim only after the specific geometry passes basis/grid/domain convergence, material-range checks, conservation or residual checks, and comparison with an independent solver, analytic result, or experiment. Visual smoothing never creates additional computed samples.

## Training-derived usability and visualization backlog

| Improvement | Status | Concrete scope |
|---|---|---|
| Phase throughout the application | Implemented within scope below | Add physically defined phase outputs wherever a complex field or scattering amplitude is available: propagating diffraction orders, reflection/transmission spectra, field maps, guided modes, resonance evidence, and exported data. Label the reference convention, phase wrapping, and locations where phase is undefined because amplitude is numerically zero. |
| TMM wavelength workflows | Implemented within scope below | Give the transfer-matrix workspace a first-class wavelength sweep with wavelength on the horizontal axis, matching the spectrum workflow and allowing phase, R, T, A, and conservation to be displayed together or selected independently. |
| Fourier-order failure recovery | Implemented within scope below | When a scan contains failed Fourier-order checks, report the failed wavelength/angle coordinates in a clickable table, distinguish exploratory points from validated points in the plot, disable figure export until important points pass, and provide an in-place checklist: open the failing single-wavelength point, raise Fourier budget, move slightly away from a grazing cutoff, and rerun. |
| Layer labels in field maps | Implemented within scope below | Draw persistent labels for each finite layer, its thickness, and actual material name on field maps. Keep interfaces legible without obscuring the data. |
| Material-specific cladding labels | Implemented within scope below | Replace generic “Cladding” labels in guided-mode field profiles with the actual upper and lower material names. State which side contains an evanescent tail when asymmetric claddings produce one. |
| Layout and scrolling | Implemented within scope below | Reduce unnecessary scrolling after the “Resonance evidence: Center, linewidth, and off-resonance fields” card through compact result summaries, collapsible detail panels, persistent operation navigation, and aligned controls. Validate with a representative desktop viewport. |
| Silicon nitride material coverage | Implemented within scope below | Ensure silicon nitride appears consistently in every relevant material dropdown, preset, project import path, and source-aware material selector. Add a regression check for this catalogue contract. |
| Material color semantics in 3D | Implemented within scope below | Assign each distinct material a stable, distinguishable color in all 3D visualizations, retain a visible material legend, and preserve contrast for color-vision accessibility. |
| Numerical-solver guidance | Implemented within scope below | Each solver should explain how to resolve its own warnings and quality gates, including convergence, Fourier orders, grazing cutoffs, boundary-domain sensitivity, and feature-fit identifiability. Link warning text directly to the relevant controls and a short worked example. |


## Training upgrade delivery — 2026-09-28

The preceding backlog descriptions retain the original requests. Delivered scope and limits:

- Persistent workspace/operation navigation and a run selector remove the need to scroll to start a calculation. Plot appearance is collapsed by default; desktop input controls scroll independently of results. The existing focus/hide controls remain available.
- A bounded, single-worker background queue serves browser calculations. Scans, band paths, adaptive rounds, Bayesian draws, bootstrap fits, and LDOS/tolerance sweeps report checkpoints. Every queued request shows elapsed time; uninstrumented matrix/optimization stages use indeterminate progress. Remaining times are estimates, and cancellation waits for the next checkpoint. Jobs are local and are not durable across server restarts.
- Scalar angle/wavelength and k-space grids support 201 points per axis; coherent polarization grids support 101; spectra support 2,001 wavelengths. Users set a smaller total scan limit in the toolbar. Bootstrap accepts up to 200 replicates. Rectangular-mode grids allow 120,000 cells; the previously failing 0.01375 µm training geometry completed on a 200 × 176 grid. This verifies execution, not mode convergence.
- Coherent TMM and specular FMM reflection/transmission phase is returned and plotted against wavelength. E/H phase maps, vector-mode phase maps, and local-basis Jones phase maps are available; zero-amplitude phase is masked. Existing Bloch phase remains available. Phase convergence is not inferred from power convergence, and spectra of incoherent power and LDOS have no single phase.
- Failed spectral coordinates have a recovery panel and single-wavelength shortcuts. Numerical figure gates remain in place. Adaptive single-feature fits now withhold Q when the linewidth/center is constrained, undersampled, or uncertain; the earlier joint-fit boundary guard is included.
- Layer labels include both background and feature materials. Planar cladding labels identify material keys. Material colors and a dynamic 3D legend replace shared fallback colors. SiN is present in the common material catalogue/selector. A contrast warning identifies patterned layers whose background and feature are identical.
- Measured/simulated overlay annotation overlap is removed and simulated lines are dashed. Parameter-sweep plots now show numerical ticks. Mode-run history records actual grid, input mesh/padding, and effective-index differences; the three-run workflow varies mesh and domain separately. It does not automatically prove mode identity.
- Normalized collected decay is exposed alongside LDOS collection. FIRST_STEPS.html is a standalone, printable introduction with a published MPB structure, analytical film control, original band illustrations, stepwise exercises, and explicit provenance. Teaching geometries and session observations are labeled separately from published reference values.

Commercial-interface references: COMSOL model/settings/graphics organization (https://www.comsol.com/support/learning-center/article/34891/12) and Lumerical visible simulation controls (https://optics.ansys.com/hc/en-us/articles/36952912384403-Ansys-Lumerical-FDTD-Modern-User-Interface). These inform interaction design; no solver equivalence or superiority claim is made.
