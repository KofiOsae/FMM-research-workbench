# FMM Research Workbench master roadmap

Status date: 2026-10-13. This roadmap consolidates the approved scientific,
device, interface, validation, tutorial, and deployment work. Detailed status
and acceptance criteria live in `TASK_TRACKER.md`; evidence lives in
`VALIDATION_MATRIX.md`.

## Completion contract

The project is complete only when every retained task is `VERIFIED`, critical
regressions pass, benchmark evidence is documented, device interfaces and
exports work, tutorials cover each supported workflow, and no critical defect
remains. Verification is always limited to the solver scope stated in the task.

## Dependency order

1. **Foundation and provenance** — versioned scientific envelopes, application
   and numerical environment manifests, resource controls, durable checkpoints.
2. **Independent component workspaces** — multilayers, Fabry–Pérot cavities,
   DBR/Bragg filters, guided modes, periodic scattering, and ports.
3. **Integrated photonics foundation** — shared complex E/H port contract,
   reference planes, mode identity, and component-to-component connections.
4. **Finite devices** — scalar grating workflow completion, then a separately
   validated full-vector propagator; no scalar-to-vector relabeling.
5. **Advanced eigensystems** — Bloch symmetry/tracking, patterned poles, QNMs,
   BIC radiation coefficients, and finite-crystal transmission.
6. **Design and uncertainty** — reusable observables, constraints, Pareto and
   robust design, fabrication rules, and quantity-specific Trust validation.
7. **Publication workflow** — common plot-data schema, comparison objects,
   citations, research bundles, and complete worked tutorials.

## Workstreams

### A. Architecture, reproducibility, and deployment

- `PM-001` establishes persistent roadmap and recovery files.
- `ARCH-001` separates solver input, workspace state, result, validation, and
  environment in a versioned project/result schema.
- `ARCH-002` exports commit, platform, Python, numerical packages, and linear
  algebra backend.
- `PLOT-001` defines one raw-data contract for scalar, complex, vector, masks,
  uncertainty, normalization, units, and reference phase.
- `DEPLOY-001` maintains bounded jobs, resource preflight, actionable failures,
  responsive health checks, and local/public limits.

### B. User interface and research workflow

- `UI-001` keeps one active tool visible by default with persistent run actions.
- `UI-002` supports desktop, laptop, tablet, and narrow screens without hiding
  scientific meaning or critical errors.
- `UI-003` applies consistent card, legend, map, hover, layer/material, and
  phase conventions.
- `UI-004` explains every parameter's meaning, units, supported range,
  numerical cost, experimental choice, and recovery action.
- `DOC-001` and `DOC-002` provide literature-based first use and a
  design → optimize → tolerance → validate → export flow for every capability.

### C. Multilayers, cavities, and filters

- `STACK-001` provides an independent multilayer stack editor with physical
  top-to-bottom ordering, material provenance, phase-aware results, and saved
  components.
- `TMM-001` provides wavelength-resolved coherent R/T/A and phase with analytic
  film controls.
- `FP-001` adds Fabry–Pérot resonance, free spectral range, finesse, linewidth,
  field buildup, loss, and mirror decomposition with stated approximations.
- `DBR-001` adds quarter-wave and arbitrary DBR/Bragg design, stop-band,
  penetration-depth/group-delay, defect-cavity, tolerance, and validation flows.
- `POLE-001` retains the uniform multilayer complex-pole solver with frozen-index
  limitations and later adds dispersive pole support only with validation.

### D. Periodic FMM/RCWA and metasurfaces

- `FMM-001` maintains order-resolved scattering, complex amplitudes, phase,
  fields, material-aware geometry, conservation, and convergence.
- `FMM-002` keeps the installed factorization scope explicit and requires
  convergence for discontinuous high-contrast media.
- `GEO-001` supports predefined and arbitrary raster/drawn/imported periodic
  structures, transforms, duplicate copies, extrusion, labels, and export.
- `META-001` builds local-periodic phase/amplitude libraries and flat-optic
  assignment with efficiency, sampling, and validity warnings.
- `OPT-001` supplies bounded multi-parameter observables and constraints;
  `OPT-002` adds Pareto reporting; `OPT-003` covers robust/yield optimization.

### E. Guided modes and integrated photonics

- `MODE-001` maintains analytic and numerical multilayer planar modes.
- `MODE-002` completes full-vector cross-section modes for imported/raster
  geometries, tensor materials, open boundaries, identity tracking, and mesh,
  padding, and mode-count convergence.
- `PORT-001` defines power-normalized complex E/H ports, direction, normal,
  reference plane, mesh, mode ID, and interpolation error.
- `IP-001` creates an integrated-photonics component workspace around those
  ports; `IP-002` connects components through modal scattering data.
- `RING-001` adds straight/bent waveguide and ring/bus eigenmode, coupling,
  resonance, FSR, Q, loss, critical-coupling, and tolerance workflows with
  literature and analytic controls.

### F. Finite grating couplers

- `FG-001` preserves the validated scalar 2D TE finite grating scope and closed
  channel budget.
- `FG-002` completes common-reference-plane phase de-embedding and an
  independently varied lateral-window validation family.
- `FG-003` adds finite top-oxide/air interfaces and matched literature geometry.
- `FG-004` connects the vector port contract to a reviewed full-vector finite
  propagator with modal reflection/transmission and radiation decomposition.
- `FG-005` addresses finite width/3D fiber overlap only through an actual 3D
  solver or an explicitly named approximation with independent evidence.

### G. Bands, poles, resonances, and emission

- `BAND-001` maintains full reciprocal-cell gap screening.
- `BAND-002` adds complex Bloch fields, phase gauge, symmetry, degeneracy
  subspaces, and crossing-aware tracking.
- `BAND-003` links an exact primitive cell to finite-crystal transmission.
- `QNM-001` adds patterned outgoing-wave poles and normalized QNMs;
  `BIC-001` derives ideal-BIC symmetry and polarization charge from isolated
  modal radiation coefficients.
- `EMIT-001` extends planar LDOS toward named-mode beta only with normalized
  modes; `EMIT-002` reserves patterned Green functions and patterned collection.

### H. Validation and publication

- `VAL-001` maintains analytic and regression controls for films, modes,
  diffraction, bands, finite devices, conservation, and convergence.
- `VAL-002` adds independent external-solver/literature cases with matched
  geometry, observable, material, polarization, and numerical conventions.
- `VAL-003` enforces quantity-specific publication thresholds and keeps
  convergence separate from benchmark agreement.
- `PUB-001` creates a complete publication bundle with inputs, outputs, units,
  versions, citations, failed points, validation, limitations, and raw figures.

## Scope rules

- Infinite periodic RCWA is never described as a finite coupler or finite
  crystal. Uniform TMM has no lateral diffraction. Bound-mode solvers do not
  establish excitation efficiency. Driven spectral fits are not poles or QNMs.
- The finite grating solver is currently scalar 2D TE and invariant across
  width. Literature headline efficiencies are trend context until the modeled
  stack, dimensionality, source, reference planes, mesh, and efficiency
  definition match.
- UI availability never implies scientific verification. Each result must carry
  its method, physical scope, assumptions, units, numerical settings, and
  validation state.

