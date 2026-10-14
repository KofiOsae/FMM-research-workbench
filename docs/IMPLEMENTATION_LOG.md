# Implementation log

## 2026-10-14 — UI-006, GEO-002, and PLOT-003

- Made the research header and workspace navigation persistent, retained all
  three working regions from 1180 px upward, and added saved Configuration and
  Results & diagnostics width controls.
- Moved the newest run status, metrics, convergence evidence, warnings, and
  export actions into the diagnostics region while keeping the active plot in
  the scientific region.
- Removed the duplicate main-field color scale, isolated line and heatmap hover
  state, and preserved the correct plot model during failed-run recovery.
- Made primitive and imported-mask previews honor solver feature offsets.
  Added scaled transformed copies, per-component raster movement, physical gap
  entry, and smoother circular preview tessellation.
- Updated the arbitrary-geometry tutorial with the exact transform and
  connected-component behavior.
- Full suite: 106 tests passed; the final HTTP asset regression passed. Live
  1265 px browser checks confirmed the three regions, adjustable widths,
  offset custom-mask preview, and a clear console.
- Committed and pushed as `040e442`. Local and Render served release marker
  `persistent-research-shell-2026-10-14`; Render also served the new geometry,
  width-control, and sticky-shell assets.

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
- Committed as `f327996`, pushed to `main`, restarted the local service, and
  verified that local and Render both serve release marker
  `three-region-workspace-2026-10-13` and the new JavaScript/CSS assets.

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

