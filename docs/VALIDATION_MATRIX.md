# Validation matrix

This matrix separates numerical regression/convergence from physical benchmark
agreement. Passing one column does not imply passing the other.

| Capability / tasks | Numerical and regression evidence | Independent physical evidence | Remaining gate |
| --- | --- | --- | --- |
| PM-001 | Protocol-file contract test: required files, IDs, states, resume fields; 106-test suite | Repository review against mandatory protocol | Preserve contract as roadmap evolves |
| ARCH-001/ARCH-002 | Version-2 sections, v1 migration, fingerprint invariance/tamper rejection, environment fields, browser save | Local browser creates a fingerprinted v2 project; environment endpoint reports runtime identity | ARCH-001 verified; Render reopen and publication-bundle reconstruction remain under ARCH-002/PUB-001 |
| UI-001/UI-002/UI-005 | HTTP assets, DOM markers, responsive CSS, run interaction | Visual desktop/laptop/tablet review of Configuration, Scientific view, Results & diagnostics | Maintain with each dynamic card |
| TMM-001/STACK-001 | Conservation, phase, project round trip | Fresnel and complex Fabry–Pérot formula | Independent stack workspace pending |
| FP-001/DBR-001 | Resonance and stop-band convergence | Analytic symmetric cavity and quarter-wave DBR | Not started |
| FMM-001/FMM-002 | Order/grid convergence, diffraction power, strict conservation | TMM uniform limit and literature/material controls | Continue case-specific convergence |
| GEO-001/META-001 | Raster orientation, copy/rotation, field/3D semantics | Local-periodic design consistency | Finite-aperture validation pending |
| MODE-001 | Root/mesh convergence and bound-mode checks | Analytic asymmetric slab | Verified within lossless 1D scope |
| MODE-002/PORT-001 | Mesh/domain tracking, self-overlap, directionality | Analytic slab and identical imported mode | General raster/PML/vector port schema |
| FG-001 | Closed budget, mesh/absorber/padding/source/interface studies | Uniform guide and analytic port controls | Verified scalar 2D TE scope |
| FG-002 | Phase de-embedding and lateral-window family | Reciprocal two-source coefficient equality | In progress |
| FG-003/FG-004 | Geometry, conservation, reciprocity, decomposition | Matched Marchetti or another declared case plus independent solver | Not started |
| BAND-001 | Basis/full-zone refinement, reciprocal frequencies | MPB square-rod reference and folded uniform waves | Verified frequencies/gap screen |
| BAND-002/BAND-003 | Gauge, subspace tracking, path reversal, finite transmission | Analytic parity and MPB symmetry | In progress/not started |
| OPT-001/OPT-003 | Feasibility, repeatable seeds, final Trust checks | Known synthetic optimum and fabrication distribution controls | Pareto frontier pending |
| EMIT-001 | Quadrature/cutoff refinement and homogeneous normalization | Wylie–Sipe/Hartman planar limits | Named-mode beta pending |
| QNM-001/BIC-001 | Pole/residue/basis and loop convergence | Independent QNM/eigenmode reference | Not started/blocked |
| VAL-001/VAL-003 | 110-test critical suite and profile regressions | Per-row benchmark evidence above | Apply gates uniformly |
| DOC-001/DOC-002 | Link/content/asset assertions | User follows expected workflow/results | Advanced worked examples pending |
| PUB-001/PLOT-001 | Export schema and figure reconstruction | Rebuilt plot equals displayed samples | Not complete |

## Critical regression command

```powershell
.\.venv\Scripts\python.exe -m unittest -v
```

For expensive new solvers, add a fast analytic control to the default suite and
document any slow benchmark command separately. Never loosen a tolerance merely
to make a regression pass without physical and numerical justification.

