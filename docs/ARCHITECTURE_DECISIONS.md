# Architecture decisions

## ADR-001 — Preserve solver boundaries

**Decision:** Periodic RCWA, uniform TMM, bound modes, uniform-port projection,
bands, scalar finite gratings, and planar Green functions remain separate
scientific workspaces with explicit scopes.

**Reason:** Similar plots do not imply compatible boundary conditions,
normalization, dimensionality, or observables. Combining labels would create
false device-level claims.

## ADR-002 — Versioned scientific envelope

**Decision:** New saved projects and publication bundles will separate
`solver_input`, `workspace_state`, `result`, `validation`, and `environment`.
Legacy imports will migrate explicitly.

**Reason:** Layout or panel changes must not alter solver fingerprints, and a
result must retain the exact physical/numerical model that created it.

## ADR-003 — No visualization-only scientific placeholders

**Decision:** A scientific control becomes user-facing only with implemented
physics, exports, documented assumptions, acceptance tests, and benchmark
evidence. Blocked 3D, QNM, patterned-emission, and BIC tasks remain absent.

## ADR-004 — Scalar finite grating remains scalar

**Decision:** The existing finite grating is named and validated as invariant-
width scalar 2D TE FDFD. Full-vector and finite-width work are separate tasks.

**Reason:** Copying scalar outputs into vector-looking plots would not provide
modal power, polarization, finite-width radiation, or 3D fiber efficiency.

## ADR-005 — Focused tool deck is the default interaction model

**Decision:** Show one active calculation card and its actions by default;
retain All tools and More for overview and advanced controls.

**Reason:** Users can edit and run without traversing a document-length page,
while complete scope, validation, and project controls remain accessible.

## ADR-006 — Convergence and benchmark agreement are separate

**Decision:** Reports and future schemas carry distinct numerical-stability and
physical-benchmark statuses.

**Reason:** A stable discretization can solve the wrong physical model, and a
single matching data point can be numerically unresolved.

## ADR-007 — Shared port contract precedes integrated-device composition

**Decision:** Integrated photonics and full-vector finite devices depend on a
serialized, power-normalized complex E/H port contract with normal, direction,
reference plane, mesh, mode ID, and interpolation error.

**Reason:** Coupling efficiency and phase cannot be transferred defensibly
between otherwise independent eigensolvers and propagators without a common
normalization and coordinate convention.

