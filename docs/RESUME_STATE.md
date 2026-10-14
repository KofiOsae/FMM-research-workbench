# Resume state

Updated: 2026-10-14

- **Current task:** `ARCH-001` — versioned project envelope and solver-input
  fingerprint separation.
- **Stage:** inspect the current project import/export and serialization paths,
  define the migration boundary, then add compatibility and fingerprint tests.
- **Most recent verified checkpoint:** `UI-006`; `GEO-002`, `PLOT-002`, and
  `PLOT-003` are tested within their stated main-preview/main-canvas scope.
- **Files changed at this checkpoint:** `app.py`, `index.html`, `workbench.js`,
  `workbench.css`, `METASURFACE_GUIDE.html`, `test_workbench_upgrade.py`, and
  this tracker state.
- **Last successful tests before this checkpoint:** 106 tests passed after the
  persistent shell, adjustable three-region interface, geometry transform,
  offset-preview, and plot-state changes. HTTP assets and live 1265 px browser
  checks passed; no browser console errors were observed.
- **Known failures:** none in tracked code. A system PowerShell profile emits a
  harmless Conda `Invoke-Expression` warning; use `.venv\Scripts\python.exe`.
- **Uncommitted changes:** inspect `git status`; do not add the user's untracked
  PDFs, PowerPoint, MATLAB, or JSON project files.
- **Immediate next action:** map the legacy version-1 project fields and the
  chained specialist-workspace collectors, then introduce a version-2 envelope
  with explicit `solver_input`, `workspace_state`, `result`, `validation`, and
  `environment` sections plus a version-1 migration test.
- **Dependencies:** preserve old JSON project import; layout/UI state must not
  change the solver-input fingerprint.
- **Decision:** build shared schemas before independent multilayer and
  integrated-photonics workspaces. Do not expose inactive full-vector, QNM,
  patterned Green-function, or 3D controls.

