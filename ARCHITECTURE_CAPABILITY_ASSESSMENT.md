# FMM Research Workbench architecture and capability gap assessment

Status: 2026-10-09. This document is the review gate before the next major
scientific solver milestone. A capability is considered implemented only when
its physical scope, normalization, validation, interface, export contract, and
regression tests are all present.

## 1. Current architecture

### Request and execution path

`workbench.js` submits numerical operations to `app.py`. Browser calculations
enter a bounded, single-worker queue in `run_jobs.py`; GET requests remain
responsive while one numerical job runs. The server dispatches each operation
to a solver-specific module and attaches the submitted input fingerprint,
software versions, method provenance, and receipt time to the response.

This architecture is appropriate for a single-user local research tool and a
small shared demonstration. It is not a distributed compute service. The
public container has lower safe memory limits than a research workstation.

### Existing solver boundaries

| Component | Existing implementation | Reusable scientific asset | Boundary that must remain explicit |
| --- | --- | --- | --- |
| Periodic scattering | `stack.py`, `grcwa` | Layer parsing, material grids, order-resolved power and complex fields, convergence studies | Infinite periodic unit cell; it is not a finite coupler or finite crystal |
| Uniform multilayers | `tmm.py` | Independent R/T/A and phase control | Laterally uniform films only |
| Planar guided modes | `waveguide.py`, `multilayer_modes.py` | Analytic asymmetric-slab roots and numerical multilayer modes | Bound 1D modes; no finite device propagation |
| Vector cross-section modes | `vector_modes.py` | Six-component fields, mode metrics, mesh/domain comparison | Isotropic built-in cores; no general imported raster or open-vector port |
| Uniform-port projection | `mode_coupling.py` | Reciprocity overlap, forward/backward mode amplitudes, source screens | One uniform plane; it is not grating coupling |
| Finite grating device | `finite_grating.py` | Normalized 2D TE port launch, finite radiation, Gaussian overlap, closed budget, spectrum/design/tolerance/validation | Scalar, lossless, uniform grating with one substrate and no finite lateral width |
| Bands | `bands.py` | TE/TM eigenfrequencies, path fields, full sampled reciprocal-cell gap screen | Closed, lossless 2D crystal; no finite-crystal transmission link yet |
| Resonance analysis | `resonance.py`, `multi_resonance.py`, `resonance_fields.py` | Adaptive driven-spectrum fits, field evidence, branch assignment | Phenomenological driven response; not a patterned pole or QNM solver |
| Optimization | `optimization.py`, `observables.py`, `research_validation.py` | Named observables, explicit constraints, robust samples, final Trust checks | Derivative-free bounded search; no adjoint/topology engine |
| Uncertainty and fitting | `tolerance.py`, `multifit.py`, `bayesian.py` | Fabrication sampling, covariance, posterior diagnostics and prediction | Conditional on entered distributions, likelihood and optical model |
| Emission | `cavity.py`, `dipole_ldos.py` | Q/V estimate and planar Green-function LDOS/collection | No patterned Green function or named patterned-mode beta factor |

## 2. Scientific integrity assessment

### Strengths already reusable

- Physical workspaces are separated: periodic FMM, uniform TMM, guided modes,
  uniform-port projection, bands, and finite grating devices do not share a
  misleading common result label.
- Solver results carry exact submitted inputs, SHA-256 fingerprint, software
  versions, method references, units in field names, and numerical settings.
- Periodic scattering reports diffraction orders, energy balance, Fourier and
  geometry-grid changes. Finite devices report every modeled power channel and
  mesh, absorber, and padding sensitivity.
- Named observables and constraints can be reused by periodic sweeps,
  optimization, tolerance studies, and final Trust validation.
- The test suite contains analytic film/mode limits, material references, an
  MPB band benchmark, diffraction conservation, modal self-overlap, tolerance,
  optimizer, and finite-device controls.

### Gaps

1. Project files still combine UI state and scientific results in one document.
   The next schema must contain separate `solver_input`, `workspace_state`,
   `result`, `validation`, and `environment` sections with a schema version.
2. Software version reporting covers numerical packages but needs application
   commit, operating system, BLAS/LAPACK provider, and solver-specific version.
3. Several plots export raw data, but there is no single plot-data contract for
   axes, units, complex values, normalization, masks, and error bars.
4. Some convergence defaults use an absolute change without a quantity-specific
   relative tolerance or uncertainty floor. Numerical convergence and physical
   benchmark agreement must remain separate report fields.
5. The local preview server cannot guarantee multi-user isolation, durable jobs,
   or recovery of a running calculation after a process restart.

## 3. Reliability findings and immediate controls

The audit found two concrete application faults and two deployment risks:

- A high-aspect-ratio lattice could make reciprocal-order index differences
  exceed the geometry FFT raster and raise an internal `IndexError`. A preflight
  now reports the minimum required geometry grid and corrective action.
- Finite-device tolerance exports contained NumPy Boolean and infinite bound
  values that strict JSON rejected. Returned types and open bounds are now
  serialized explicitly.
- Completed queue results could retain large field arrays for 30 minutes. The
  browser now acknowledges consumed results so the server releases them while
  durable summaries remain.
- Direct synchronous requests could run beside the queued solver. A process-wide
  solver gate now prevents concurrent large numerical jobs; the HTTP listen
  backlog is enlarged so health and polling requests remain responsive.

For the shared Render demonstration, numerical libraries are restricted to one
thread and the finite FDFD grid is capped at 35,000 cells. Local installations
retain the 120,000-cell ceiling. These controls reduce process termination and
503 responses; they do not turn the free shared service into a high-performance
compute platform.

## 4. Capability gaps against the requested research chain

### Guided mode → finite device → far field

The scalar finite grating already uses its own normalized TE port eigenfield and
decomposes residual forward mode, backward guided reflection, upward radiation,
substrate radiation, Gaussian target overlap, and power residual. The missing
general chain is a shared **vector port contract** accepted by both the vector
eigensolver and a full-vector finite-domain propagator. It must preserve complex
E/H fields, port normal, propagation direction, normalization power, mode ID,
and mesh interpolation error.

The validated silicon example currently covers a uniform scalar 2D TE device.
A quantitative reproduction of the 260 nm SOI result from Marchetti et al.
(2017) additionally requires per-tooth period/fill apodization, finite oxide,
silicon handle, and matched full-vector physics. The reported 83% value must not
be used as an acceptance target for a different geometry.

### Band eigenproblem → finite crystal

Complete sampled Brillouin-zone gap screening exists. Complex Bloch field export
is partial, while symmetry labels, degeneracy subspaces, gauge-consistent mode
tracking, and a finite-crystal transmission builder are missing. A reliable link
must create a finite crystal from the exact primitive cell and compare attenuation
inside a complete gap without identifying a symmetry-path gap as complete.

### Analysis and comparison

Complex diffraction amplitudes, phase, signed fields, vertical Poynting flux,
modal overlap, angular maps, linked sweeps, and raw-data exports exist in several
workspaces. They need a common scientific plot contract and a comparison object
that records interpolation, normalization, reference phase, and incompatible
model warnings.

## 5. Staged implementation plan and acceptance tests

### Milestone 0 — reliability and reproducibility foundation

Deliverables:

- Versioned scientific result envelope separated from workspace metadata.
- Application commit and numerical backend details in every export.
- Queue result acknowledgement, resource preflight, public-service limits, and
  actionable recovery messages.
- Capability matrix and solver-selection language in the tutorial.

Acceptance tests:

1. Every API result round-trips through strict JSON with `allow_nan=False`.
2. A project can change panel layout without changing the solver-input hash.
3. Excessive grid/order/domain requests fail before allocation and identify the
   controlling inputs and safe corrective action.
4. Static pages and `/health` remain responsive during one queued calculation.
5. Repeated completed field jobs do not leave their full arrays in server memory
   after the browser receives them.

### Milestone 1 — finite grating research workflow

Deliverables:

- Illustrated geometry with coordinate system, tooth/trench definition,
  material regions, port, radiation monitors, Gaussian target, angle and waist.
- Contextual definitions and valid ranges for every parameter.
- Result interpretation panel with energy accounting, good/review/fail gates,
  dominant loss channel, and specific design actions.
- Per-tooth uniform, linear-apodized, and imported grating tables; finite oxide
  and substrate layers.
- Shared normalized TE/TM vector-port data model, then full-vector propagation
  as a separate reviewed milestone.
- A matched literature case and independent solver comparison.

Acceptance tests:

1. Uniform waveguide returns at least 97% residual forward guided power and less
   than 3% budget residual at the documented screening mesh.
2. Analytic and numerical port indices converge toward one another under mesh
   refinement; the comparison is not hidden behind a loose pass label.
3. For every solved device, modeled channel power plus numerical residual equals
   incident power within the selected tolerance.
4. Reciprocity gives the same selected-mode efficiency in both directions for a
   reciprocal matched test.
5. Target overlap never exceeds total upward radiation beyond tolerance.
6. A literature reproduction uses the same layer stack, tooth table, wavelength,
   polarization, angle, fiber convention and efficiency definition; otherwise it
   is labeled a trend comparison.
7. Spectrum, optimizer, tolerance and Trust reports preserve the same device and
   port fingerprints.

### Milestone 2 — band symmetry and mode tracking

Deliverables:

- Complex Bloch E/H field export with a documented phase gauge.
- Little-group symmetry operations and parity/character estimates at compatible
  high-symmetry points.
- Degeneracy clustering and subspace-aware tracking through crossings.
- Controlled lattice-symmetry-breaking experiment templates.

Acceptance tests:

1. MPB reference frequencies remain within the existing tolerance.
2. Symmetry eigenvalues agree with analytic parity for uniform and simple rod
   controls.
3. Rotations within a degenerate subspace do not change the reported subspace ID.
4. Forward and reverse k paths return reciprocal frequencies for reciprocal media.
5. Tracking is stable under path refinement and reports ambiguity rather than
   forcing a branch assignment.

### Milestone 3 — Bloch-to-finite-crystal transmission

Deliverables:

- Finite repetition builder from the exact band primitive cell.
- Compatible source/port selection and transmission/reflection spectrum.
- Comparison view linking Bloch bands, projected k, finite transmission and
  field decay.

Acceptance tests:

1. A homogeneous cell gives the analytic folded dispersion and no artificial
   stop band in finite transmission.
2. Increasing finite periods deepens attenuation inside a validated complete gap.
3. Interface termination is explicit and changing it is recorded as a physical
   model change.
4. Energy conservation and basis/grid convergence pass independently.

### Milestone 4 — unified scientific comparison and publication package

Deliverables:

- Common raw plot-data schema for scalar, complex, vector and uncertainty data.
- Side-by-side comparison with explicit alignment and normalization.
- Publication bundle containing inputs, outputs, provenance, failed points,
  convergence, benchmark status, limitations and figure data.

Acceptance tests:

1. Every visible quantitative plot has downloadable numeric axes and values.
2. Comparing incompatible normalization or boundary conditions produces a clear
   warning and cannot silently compute a difference.
3. Rebuilding a figure from the exported raw data reproduces plotted values.
4. Numerical convergence and benchmark agreement have separate status fields.

## 6. Delivery rule

Each milestone is independently reviewable and keeps previous project imports
working through explicit schema migration. New controls first expose the physical
definition, then the solver, then validation, and finally optimization. A plot or
UI control alone does not count as a scientific capability.

Reference implementations: MPB documentation for band eigenproblems and symmetry,
and Meep documentation for eigenmode sources, mode decomposition, diffraction
orders, flux normalization, and convergence studies. Their presence is design
guidance, not evidence that this Workbench produces equivalent results.
