# FMM Research Workbench development protocol

These repository instructions apply to every Codex session.

## Start of a session

1. Inspect `git status` and the latest commits.
2. Read `docs/MASTER_ROADMAP.md`, `docs/TASK_TRACKER.md`,
   `docs/IMPLEMENTATION_LOG.md`, `docs/RESUME_STATE.md`,
   `docs/VALIDATION_MATRIX.md`, and `docs/ARCHITECTURE_DECISIONS.md`.
3. Verify the checkpoint against the repository. Resume an interrupted task or
   select the next unblocked task by dependency and priority.
4. Preserve independently developed and untracked user files.

## Delivery rules

- Use stable task identifiers and only these states: `NOT_STARTED`,
  `IN_PROGRESS`, `BLOCKED`, `IMPLEMENTED`, `TESTED`, `VERIFIED`.
- Mark a task `VERIFIED` only after its acceptance criteria, tests,
  documentation, and supported-scope statement are complete.
- Keep periodic, uniform multilayer, guided-mode, band, finite-device, and
  emission solver boundaries explicit. Never relabel an approximation as a
  rigorous full-vector or finite-domain solver.
- Implement the smallest complete scientific capability: physical definition,
  solver, interface, export, validation, documentation, and regression test.
  A plot or disabled placeholder does not count.
- Separate numerical convergence from agreement with an analytical,
  independent-solver, literature, or experimental benchmark.
- Preserve backward compatibility or add an explicit schema migration.
- Run focused tests while developing and the critical regression suite before
  marking work verified.
- Update `TASK_TRACKER.md`, `IMPLEMENTATION_LOG.md`, `RESUME_STATE.md`, and
  `VALIDATION_MATRIX.md` at meaningful checkpoints. Record blockers and move to
  an independent task when possible.
- Commit completed reviewable checkpoints when repository workflow permits.

## End of a session

Record the current task, stage, changed files, last successful tests, failures,
uncommitted changes, immediate next action, dependencies, and decisions in
`docs/RESUME_STATE.md`. Report completed, active, blocked, and next tasks
without describing unfinished work as complete.

