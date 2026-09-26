# FMM Research Workbench v0.1.0 — research preview

FMM Research Workbench is an open-source browser interface for periodic and layered photonics calculations. It combines Fourier modal scattering, transfer-matrix checks, guided modes, two-dimensional photonic bands, resonance analysis, measured-data fitting, parameter sweeps, optimization, planar dipole LDOS, material management, visualization, and guided tutorials.

## Included workflows

- Finite-stack RCWA/FMM spectra, diffraction orders, fields, angle–wavelength maps, and incident k-space maps
- Uniform-stack transfer-matrix calculations and FMM comparison
- Planar and full-vector waveguide modes
- Infinite two-dimensional photonic-crystal bands and scalar eigenfields
- Adaptive single- and multi-resonance fitting, uncertainty and branch tracking
- GMR, quasi-BIC, plasmon, far-field polarization, and winding diagnostics with documented limits
- Parameter sweeps, constrained optimization, measurement fitting, tolerance studies, and Bayesian inference
- Planar dipole LDOS and objective collection calculations
- Dispersive material catalog and local measured-material import
- Project export, publication figures, methods reports, tutorials, and literature references

## Verification status

The release includes 56 automated physics and regression tests. These cover analytic limits, energy balance, TMM/FMM agreement, material dispersion, field geometry, band benchmarks, guided modes, resonance analysis, uncertainty workflows, and planar LDOS normalization.

Passing tests and numerical convergence do not validate every user-defined structure. Research claims require case-specific convergence, correct material provenance, and comparison with suitable independent evidence. See `VALIDATION.md` and `RESEARCH_READINESS.md`.

## Public demonstration limits

The free demonstration is intended for evaluation and tutorials. Compute-intensive jobs may be slow or limited. Server-side material import is disabled on the shared demonstration; download the project for local material imports and larger studies.

## AI assistance acknowledgement

Development of FMM Research Workbench was substantially assisted by **ChatGPT Codex by OpenAI**, including source analysis, implementation, interface design, documentation, literature-guided test planning, and automated test development. Scientific interpretation, licensing review, validation decisions, and published results remain the responsibility of the project maintainers and users.

## License

GPL-3.0-or-later. The project uses `grcwa`, whose GPL notice and attribution are preserved in `LICENSE`.
