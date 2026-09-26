# Research readiness and ten-item roadmap

Status date: 2026-09-26. A feature is marked **complete** only when its solver, interface, tutorial, and automated reference checks cover the stated scope.

| # | Required improvement | Current status | Exact scope and remaining gate |
|---|---|---|---|
| 1 | Patterned-stack complex S-matrix poles | Missing | The available pole tool is restricted to uniform multilayers with frozen optical indices. A patterned outgoing-wave complex-frequency solver is still required. |
| 2 | Quasinormal-mode fields and normalization | Missing | No normalized patterned QNM field is returned. The field workspaces are driven-frequency fields or closed/bound eigenfields. |
| 3 | Radiative, absorptive, and total Q for patterned structures | Partial | Far-field linewidth fitting and uniform-stack loss-removal estimates exist. Patterned pole-based Q decomposition is not implemented. |
| 4 | Ideal-BIC identification and symmetry classification | Partial | Signed-angle, symmetry, convergence, Q-versus-angle/asymmetry, and field diagnostics are available. They do not rigorously certify an ideal BIC. |
| 5 | Polarization-vortex winding and topological charge | Partial | Closed-loop winding with complex sideband background subtraction is available as a stability diagnostic. Certification still requires radiation coefficients from an isolated patterned pole or eigenmode. |
| 6 | Full-vector 2D waveguide modes | Partial | A rectangular, isotropic dielectric cross section returns all six complex components using a vector finite-difference eigenproblem. Arbitrary geometry, tensor media, open/PML boundaries, and automatic mesh convergence are still required. |
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
| Tutorial covers every available workspace | Partial | The tutorial covers the main solver distinctions and new advanced cards. It still needs a searchable parameter glossary, complete worked examples for every advanced card, and synchronized acceptance tests for all instructions. |
| Every 2D quantitative map has axes, units, title, numeric colorbar, hover values, and working palette choice | Partial | The main scattering and several field maps meet most of this contract. Vector modes and several secondary maps do not yet share the complete renderer; some palette selectors recolor only the legend rather than the raster. |
| Professional task organization | Partial | Workspaces, calculation selectors, project save/open, source-aware materials, methods export, reference cases, and warnings are present. Long cards, inconsistent result controls, lack of searchable help, and lack of a background-job manager still impede routine use. |
| Every parameter explains meaning, units, valid range, numerical impact, and experimental choice | Partial | Major FMM inputs are explained. New mode, Bayesian, pole, constitutive, and optimization parameters need contextual definitions, recommended ranges, and consequence warnings beside each control. |
| Literature and method provenance | Partial | Core FMM, materials, bands, GMR, energy, BIC, vector-mode, planar Green-function, Bayesian-diagnostic, and QNM sources are listed. Each solver result still needs a direct method citation in its export manifest and additional primary references for coupled thermal, nonlinear, gain, and magnetic implementations before those solvers are built. |

## Release rule

Do not describe the application as generally research-grade or as a Lumerical replacement. A result becomes suitable for a research claim only after the specific geometry passes basis/grid/domain convergence, material-range checks, conservation or residual checks, and comparison with an independent solver, analytic result, or experiment. Visual smoothing never creates additional computed samples.
