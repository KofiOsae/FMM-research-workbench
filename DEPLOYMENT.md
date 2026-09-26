# Public release readiness

The workbench can be published online after a release and security pass. The current `app.py` server is a local research preview bound to `127.0.0.1`; it is not a production web server.

## Required before a public launch

1. Choose and add a project license. The installed `grcwa` dependency reports a GPL license, so distribution and hosted use must be reviewed for compatibility and source obligations. Preserve citations and licenses for material datasets.
2. Put the Python application behind a maintained production server and HTTPS reverse proxy. Add request timeouts, per-user compute quotas, job cancellation, process isolation, and monitoring. RCWA and band calculations can consume substantial CPU and memory.
3. Replace the shared `user_materials.json` store with per-user or per-project storage. Validate authentication and authorization if saved projects are added. Define retention and privacy rules for uploaded optical constants and experimental spectra.
4. Add CSRF protection, restrictive security headers, structured request logging without uploaded data, dependency scanning, and rate limiting. Keep the request-size and numeric bounds already enforced by the solver.
5. Run the complete validation suite in CI on every supported platform. Add browser tests for all tabs and public deployment smoke tests. Pin dependencies and publish exact version and dataset provenance.
6. Label the software as numerical research software, document model limits, and avoid presenting convergence as experimental validation. Provide a reproducible issue template with exported JSON settings.

## Release gate

A public release should require: all automated checks pass; every bundled benchmark reproduces its stated tolerance; lossy and lossless reference cases pass; uploaded data remain isolated; no unresolved high severity dependency findings; and an independent researcher can reproduce the tutorial results from a clean installation.
