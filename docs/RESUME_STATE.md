# Resume state

Updated: 2026-10-13

- **Current task:** `ARCH-001` — versioned project envelope and solver-input
  fingerprint separation.
- **Stage:** inspect the current project import/export and serialization paths,
  define the migration boundary, then add compatibility and fingerprint tests.
- **Most recent verified checkpoint:** `PM-001`, `UI-001`, `UI-002`, and
  `UI-005`; `PLOT-002` is tested within the main spectrum/field canvas scope.
- **Files changed at this checkpoint:** `AGENTS.md`, `README.md`,
  `test_workbench_upgrade.py`, and all files under `docs/` listed by PM-001.
- **Last successful tests before this checkpoint:** 106 tests passed after the
  persistent protocol, three-region interface, responsive layout, and main
  spectrum/field canvas changes. HTTP asset integration and live run checks
  passed; no browser console errors were observed.
- **Known failures:** none in tracked code. A system PowerShell profile emits a
  harmless Conda `Invoke-Expression` warning; use `.venv\Scripts\python.exe`.
- **Uncommitted changes:** inspect `git status`; do not add the user's untracked
  PDFs, PowerPoint, MATLAB, or JSON project files.
- **Immediate next action:** commit and deploy PM-001/UI-005/PLOT-002, verify
  local and Render markers, then inspect project serialization in `index.html`,
  `workbench.js`, and `app.py` for ARCH-001.
- **Dependencies:** preserve old JSON project import; layout/UI state must not
  change the solver-input fingerprint.
- **Decision:** build shared schemas before independent multilayer and
  integrated-photonics workspaces. Do not expose inactive full-vector, QNM,
  patterned Green-function, or 3D controls.

