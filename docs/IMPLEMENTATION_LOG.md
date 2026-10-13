# Implementation log

## 2026-10-13 — UI-005 professional three-region workspace

- Replaced the document-style focused page with independent Configuration,
  Scientific view, and Results & diagnostics regions on desktop.
- Kept geometry/source controls in a full-height scrollable column, added a
  scientific-view selector, kept diagnostics beside the active plot, and
  reduced the top command bar to the primary run path plus optional controls.
- Main horizontal and vertical field plots now remove stale line-chart hover
  metadata, report nearest sampled field values, and show both physical axes;
  the horizontal magnitude plot includes a numeric color scale.
- Verified the desktop three-column layout, laptop two-column layout, portrait
  tablet stack, scientific-view switching, a completed spectrum run, and a
  completed field run. The newest output remains visible without searching the
  page, and the browser console remained clear.
- Full suite: 106 tests passed. `UI-005` is verified. `PLOT-002` is tested for
  the main spectrum/field canvas; the separate audit of every secondary plot
  remains tracked under `UI-003` and `PLOT-001`.

## 2026-10-13 — PM-001 persistent execution framework

- Added repository-level instructions and six persistent roadmap/checkpoint
  documents with stable task identifiers and controlled statuses.
- Consolidated prior interface, solver, validation, finite-device, band,
  multilayer, integrated-photonics, resonance, emission, tutorial, and
  publication requests into one dependency-ordered tracker.
- Reconciled claims with `ARCHITECTURE_CAPABILITY_ASSESSMENT.md`,
  `RESEARCH_READINESS.md`, `LITERATURE_IMPLEMENTATION_MATRIX.md`, existing
  tests, and current code. Partial or missing solvers remain explicitly open.
- Added a regression contract for required files, task states, identifiers, and
  resume fields.

## 2026-10-13 — UI-001/UI-002 focused workspace

- Commit `3c393a4` added the one-card Focused deck, previous/next/All/More
  navigation, compact chrome, and persistent current-tool actions.
- Corrected portrait-tablet sticky offset and visually checked 1440×1000,
  1024×768, and 768×1024 viewports.
- Full suite: 105 tests passed. Local and Render served
  `focused-tool-deck-2026-10-10` and the new assets.

## Earlier verified foundations

- Periodic FMM/TMM controls, guided modes, full-zone bands, constrained design,
  uncertainty, planar LDOS, custom periodic geometry, finite scalar grating,
  provenance, queue/resource controls, tutorials, and deployment evidence are
  summarized in `RESEARCH_READINESS.md` and
  `ARCHITECTURE_CAPABILITY_ASSESSMENT.md`.
- Their exact scope is preserved in `TASK_TRACKER.md`; no historical “complete”
  label expands a solver beyond its tested physics.

