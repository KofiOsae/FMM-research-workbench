# Fourier-factorization audit

Audit date: 2026-10-09.

The periodic solver delegates patterned-layer Fourier algebra to the installed
`grcwa` package. Inspection of `grcwa/fft_funs.py::Epsilon_fft` shows that it
constructs a convolution matrix of the sampled permittivity and obtains an
inverse-permittivity operator by matrix inversion. `grcwa/rcwa.py` passes this
operator into `MakeKPMatrix` for patterned layers.

This is useful inverse-rule behavior, but the inspected package does not expose
an interface-normal vector field or a normal-vector factorization path. The
Workbench therefore records the implementation as **epsilon convolution plus
matrix inverse**, not as a complete Li normal-vector implementation.

The relevant mathematical reference is L. Li, “Use of Fourier series in the
analysis of discontinuous periodic structures,” *JOSA A* 13, 1870–1876 (1996),
[doi:10.1364/JOSAA.13.001870](https://doi.org/10.1364/JOSAA.13.001870).

## Consequences for results

- Fourier-order convergence is required for every quantitative patterned-layer
  result.
- TM and high-contrast discontinuous geometries may converge slowly.
- A stable low-order value is not proof of convergence.
- Geometry-grid convergence is independent of Fourier-basis convergence.
- The UI and exports must not describe the current solver as using a complete
  Li normal-vector factorization.

## Acceptance gate for a future factorization upgrade

An upgrade must document the exact direct/inverse product rules, interface
normal construction, field reconstruction, and polarization conventions. It
must reproduce analytic lamellar controls, improve convergence on discontinuous
TM benchmarks, preserve energy balance, and retain the current implementation
as a regression option until the comparison is reviewed.
