# Master task tracker

Allowed states: `NOT_STARTED`, `IN_PROGRESS`, `BLOCKED`, `IMPLEMENTED`,
`TESTED`, `VERIFIED`. Verification is restricted to the stated scope.

| ID | Status | Depends on | Acceptance criteria / evidence gate |
| --- | --- | --- | --- |
| PM-001 | VERIFIED | — | All protocol files exist, are version controlled, cross-linked, and regression-tested for required fields and valid states. |
| ARCH-001 | VERIFIED | PM-001 | Versioned project schema separates `solver_input`, `workspace_state`, `result`, `validation`, and `environment`; v1 projects migrate; fingerprint tampering is rejected; layout changes do not alter solver identity; browser save and full regression suite pass. |
| ARCH-002 | IMPLEMENTED | ARCH-001 | Exports include application commit, OS, Python, package and solver versions, and BLAS/LAPACK provider; round-trip strict JSON test passes. |
| PLOT-001 | NOT_STARTED | ARCH-001 | One documented raw plot-data schema reconstructs every visible quantitative plot and warns on incompatible comparison. |
| DEPLOY-001 | VERIFIED | PM-001 | Local and Render health/assets respond; queue/resource preflight and result release regressions pass; service limits are exported. |
| UI-001 | VERIFIED | PM-001 | Focused one-card deck, previous/next/All/More controls, and persistent active actions work; HTTP asset regression passes. |
| UI-002 | VERIFIED | UI-001 | Desktop 1440×1000, laptop 1024×768, and portrait tablet 768×1024 show usable controls without covered card headings. |
| UI-003 | IN_PROGRESS | PLOT-001 | All dynamic cards, multi-series legends, 2D maps, material colors, phase masks, hovers, and exports share the documented renderer contract. |
| UI-004 | IN_PROGRESS | UI-001 | Every advanced control explains meaning, units, valid/recommended range, cost, physical choice, and warning recovery. |
| UI-005 | VERIFIED | UI-001, UI-002 | Desktop uses Configuration, Scientific view, and Results & diagnostics regions; newest output is visible beside its quality evidence; top controls are compact; desktop/laptop/tablet visual and run interaction checks pass. |
| UI-006 | VERIFIED | UI-005 | Header/workspace navigation remains visible, three columns persist from 1180 px, configuration and diagnostics widths are user adjustable and saved, and live 1265 px visual/console checks pass. |
| PLOT-002 | TESTED | UI-005 | Main spectrum/field canvas clears stale plot metadata, field maps have both axes, numeric scales and nearest-sample hover; field execution and asset regressions pass; full secondary-plot audit remains under UI-003/PLOT-001. |
| PLOT-003 | TESTED | PLOT-002 | Main heat and signed-field maps use one color scale, clear line-chart state, preserve the correct hover model across success/failure recovery, and pass asset/syntax regressions; secondary canvases remain under UI-003. |
| STACK-001 | NOT_STARTED | ARCH-001 | Independent multilayer editor supports top-to-bottom layers, provenance, reusable components, phase-aware spectra, project migration, tutorial, and tests. |
| TMM-001 | VERIFIED | — | Coherent wavelength-axis R/T/A/phase and conservation match analytic Fresnel/Fabry–Pérot controls. |
| FP-001 | NOT_STARTED | STACK-001 | Fabry–Pérot card reports resonances, FSR, finesse, linewidth/Q, buildup, loss and mirror decomposition against analytic controls. |
| DBR-001 | NOT_STARTED | STACK-001 | Quarter-wave/arbitrary DBR and defect filters expose stop-band, delay/penetration, tolerance, optimization, and analytic/literature validation. |
| POLE-001 | TESTED | STACK-001 | Uniform frozen-index poles and Q estimates pass regressions; dispersive extension remains outside verified scope. |
| FMM-001 | VERIFIED | — | Periodic order-resolved power, complex amplitudes, phase, fields, conservation, failure recovery, and order/grid convergence tests pass. |
| FMM-002 | VERIFIED | FMM-001 | Export and tutorial state installed factorization limits; no full Li normal-vector claim; audit regression passes. |
| GEO-001 | TESTED | FMM-001 | Draw/import raster, straight lines, transforms, duplicate/rotate copy, extrusion, 3D labels/colors, and custom-mask regressions pass; fabrication-rule editor remains. |
| GEO-002 | TESTED | GEO-001 | Solver feature offsets appear in primitive and imported-mask previews; transformed copies support anisotropic scale; disconnected raster features can move independently or receive a physical horizontal gap; touching regions intentionally become one feature. |
| META-001 | IMPLEMENTED | GEO-001, OPT-001 | Local-periodic library and flat-lens guidance exist; end-to-end finite aperture/full-device verification remains. |
| OPT-001 | VERIFIED | FMM-001 | Named observables, one-to-five bounded variables, independent constraints, feasibility, final Trust validation, and regressions pass. |
| OPT-002 | NOT_STARTED | OPT-001 | Pareto frontier preserves candidates and reports efficiency/balance tradeoffs without requiring a hidden penalty weight. |
| OPT-003 | TESTED | OPT-001 | Fixed-sample robust objective, correlated tolerance, yield summaries, seeds, and regression tests pass; optimizer-wide fabrication rules remain. |
| MODE-001 | VERIFIED | — | Analytic slab and numerical multilayer modes agree under documented lossless bound-mode controls. |
| MODE-002 | TESTED | MODE-001 | Six-component built-in-core solver and mesh/domain history pass; imported raster, tensors, PML/open boundaries, and mode-count gate remain. |
| PORT-001 | TESTED | MODE-002 | Reciprocity overlap, forward/backward amplitudes, identical-mode control, and half-mesh tracking pass; shared serialized port contract remains. |
| IP-001 | NOT_STARTED | STACK-001, PORT-001 | Integrated-photonics workspace defines components, ports, coordinate/reference planes, and compatible solver selection without placeholder physics. |
| IP-002 | NOT_STARTED | IP-001 | Components connect through normalized modal scattering data with conservation, reciprocity, phase, and cascade tests. |
| RING-001 | NOT_STARTED | IP-001, MODE-002 | Ring/bus design reports modes, coupling, resonance, FSR, intrinsic/coupling/loaded Q, loss, tolerance, and analytic/literature checks. |
| FG-001 | VERIFIED | MODE-001 | Scalar 2D TE finite grating reports calibrated bidirectional modal/radiation channels, closed budget, spectrum/design/tolerance, fields, and validation. |
| FG-002 | IN_PROGRESS | FG-001, PORT-001 | De-embed complex coefficients to common reference planes and add independent lateral-window validation; phase reciprocity acceptance passes. |
| FG-003 | NOT_STARTED | FG-002 | Finite top oxide/air geometry and matched source/monitor conventions reproduce a declared literature case within stated scalar limits. |
| FG-004 | NOT_STARTED | PORT-001, FG-002 | Full-vector finite propagator accepts serialized ports and passes uniform, conservation, reciprocity, mesh, absorber, and mode-decomposition tests. |
| FG-005 | BLOCKED | FG-004 | Finite-width/3D efficiency requires a selected 3D Maxwell backend and independent resource/benchmark plan. |
| BAND-001 | VERIFIED | — | Full reciprocal-cell screening, path dispersion, fields, and MPB frequency regression pass. |
| BAND-002 | IN_PROGRESS | BAND-001, PLOT-001 | Complex field gauge, symmetry characters, degeneracy subspaces, reciprocal path and crossing-aware tracking tests pass. |
| BAND-003 | NOT_STARTED | BAND-002, IP-001 | Exact primitive-cell finite repetition shows no false homogeneous gap and increasing attenuation inside a verified complete gap. |
| QNM-001 | NOT_STARTED | FMM-002, ARCH-001 | Patterned outgoing complex poles, residues, normalized QNMs and radiative/absorptive Q pass independent controls. |
| BIC-001 | BLOCKED | QNM-001, BAND-002 | Ideal-BIC symmetry and winding use isolated modal radiation coefficients and converge under loop/basis refinement. |
| EMIT-001 | TESTED | MODE-001 | Planar LDOS/objective collection regressions pass; named-mode beta requires a normalized compatible mode projection. |
| EMIT-002 | BLOCKED | QNM-001 | Patterned Green function, patterned LDOS/beta/collection require the patterned modal/finite-domain backend. |
| VAL-001 | VERIFIED | — | Critical suite covers films, materials, modes, diffraction, bands, finite devices, conservation, optimization, uncertainty, queue and assets. |
| VAL-002 | IN_PROGRESS | VAL-001 | Every major solver has at least one matched independent analytical, solver, literature, or experimental case with explicit observable conventions. |
| VAL-003 | IMPLEMENTED | VAL-001 | Exploratory/research/publication profiles and quantity-specific gates exist; remaining work is applying them uniformly to every solver. |
| DOC-001 | VERIFIED | VAL-001 | Printable first-use guide contains literature structure, expected results, interpretation, failures, sources, and solver distinctions. |
| DOC-002 | IN_PROGRESS | UI-004 | Every capability has a synchronized design → optimize → tolerance → validate → export flow and acceptance-tested worked example. |
| PUB-001 | IN_PROGRESS | ARCH-001, PLOT-001, VAL-002 | Methods report exists; complete publication bundle and figure reconstruction test remain. |

## Current dependency frontier

The next unblocked priority is `ARCH-001`, followed by `STACK-001` and
`IP-001`. `FG-002` and `BAND-002` can proceed independently after their shared
contracts are defined. `FG-005`, `BIC-001`, and `EMIT-002` are intentionally
blocked rather than represented by visualization-only controls.

