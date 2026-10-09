# Literature-to-implementation matrix

Status: 2026-10-09. This is an auditable map from primary literature to the
actual implementation. A citation describes the governing method; it is not by
itself validation of this software.

| Capability | Primary source and role | Code path | Automated evidence | Current status and limitation |
| --- | --- | --- | --- | --- |
| Periodic Fourier modal scattering | L. Li, “Use of Fourier series in the analysis of discontinuous periodic structures,” *JOSA A* 13, 1870–1876 (1996), [doi:10.1364/JOSAA.13.001870](https://doi.org/10.1364/JOSAA.13.001870) | `stack.py`; installed `grcwa/fft_funs.py` and `grcwa/rcwa.py` | order convergence, power balance, diffraction-order regressions in `test_validation.py` | `grcwa` builds epsilon convolution matrices and uses an inverse matrix for the patterned-layer inverse-permittivity operator. The installed code does not expose an interface-normal vector-field factorization. Therefore the Workbench must not claim a complete Li normal-vector implementation; Fourier convergence remains mandatory. |
| Scalar interface averaging | A. Farjadpour *et al.*, “Improving accuracy by subpixel smoothing in the finite-difference time domain,” *Optics Letters* 31, 2972–2974 (2006), [doi:10.1364/OL.31.002972](https://doi.org/10.1364/OL.31.002972) | `finite_grating.py::_permittivity` | finite BOX/material assignment, fractional-cell, and legacy-raster tests | The scalar TE mass coefficient is cell averaged by deterministic supersampling. This is a limited scalar correction. It is not the paper’s anisotropic full-vector tensor smoothing, and the legacy binary path remains selectable for regression. |
| Finite-domain Maxwell context and validation practice | A. F. Oskooi *et al.*, “MEEP: A flexible free-software package for electromagnetic simulations by the FDTD method,” *Computer Physics Communications* 181, 687–702 (2010), [doi:10.1016/j.cpc.2009.11.008](https://doi.org/10.1016/j.cpc.2009.11.008) | `finite_grating.py`, `run_jobs.py` | mesh, absorber, padding, interface-raster, and closed-power tests | The Workbench uses its own frequency-domain scalar operator, not Meep. The paper guides source/boundary/convergence controls only. |
| SOI grating-coupler physics | D. Taillaert, P. Bienstman, and R. Baets, “Compact efficient broadband grating coupler for silicon-on-insulator waveguides,” *Optics Letters* 29, 2749–2751 (2004), [PubMed record](https://pubmed.ncbi.nlm.nih.gov/15605493/) and Taillaert *et al.*, “Grating couplers for coupling between optical fibers and nanophotonic waveguides,” *JJAP* 45, 6071 (2006), [doi:10.1143/JJAP.45.6071](https://doi.org/10.1143/JJAP.45.6071) | finite-device workspace and `finite_grating.py` | spectrum, design, tolerance, port-mode, and power-channel tests | The shared physical model now supports upper cladding, device silicon, finite BOX, and handle. It remains a uniform, invariant-width 2D scalar TE model; it cannot reproduce a complete 3D fiber-coupler efficiency. |
| Waveguide modes and modal power | A. W. Snyder and J. D. Love, *Optical Waveguide Theory* (1983); S. Fallahkhair, K. S. Li, and T. E. Murphy, “Vector finite difference modesolver for anisotropic dielectric waveguides,” *JLT* 26, 1423–1431 (2008), [doi:10.1109/JLT.2008.923643](https://doi.org/10.1109/JLT.2008.923643) | `waveguide.py`, `vector_modes.py`, `mode_coupling.py`, `finite_grating.py::_te_mode` | analytic slab and numerical-mode comparisons; overlap self-tests | The finite grating uses a scalar TE port normalized by integral |Ey|² dz and calibrates power with a reference solve. General vector E/H power normalization exists in the separate uniform-port workspace and is not yet connected to a full-vector finite propagator. |
| Reciprocity | Lorentz reciprocity as treated in Snyder and Love and standard electromagnetic scattering theory | `mode_coupling.py`; `finite_grating.py`; `finite_grating_bidirectional.py` | uniform-port overlap self-tests; two-source finite-grating coefficient comparison | The finite-grating certificate runs independent guided-mode and Gaussian equivalent-current solves and compares power-normalized coefficient magnitude and efficiency. Raw phase is exported; phase reciprocity is not certified until port reference planes are de-embedded. |
| Gaussian equivalent-current source | Total-field/scattered-field and equivalent-current formulations in Taflove and Hagness; Gaussian mode context in Taillaert et al. | `finite_grating_bidirectional.py::solve_fiber_incident` | homogeneous reference reproduction, sampled-angle report, paired-source reciprocity regression | The incident field is a 2D width-invariant Gaussian sheet synthesized from propagating angular-spectrum components and truncated inside the absorber. It is not a circular full-vector 3D fiber field. |
| Literature device target | R. Marchetti *et al.*, “High-efficiency grating-couplers: demonstration of a new design strategy,” *Scientific Reports* 7, 16670 (2017), [doi:10.1038/s41598-017-16505-z](https://doi.org/10.1038/s41598-017-16505-z) | tutorial research workflow | tutorial checks only | This reference motivates geometry, efficiency definitions, and validation. Its reported efficiency is not an acceptance value until the same tooth table, layer stack, materials, polarization, angle, and 3D model are reproduced. |

## Equations implemented in the current finite-device milestone

The time convention is exp(-i omega t). For an invariant, isotropic,
non-magnetic structure and E = y-hat Ey, the solver discretizes

`(d_x^2 + d_z^2 + k0^2 epsilon_r) Ey = b`.

The reconstructed relative magnetic field is

`Hx = i d_z(Ey)/k0`, `Hz = -i d_x(Ey)/k0`,

and the relative time-averaged flux components are proportional to

`Sx = Re(Ey conj(Hz))/2`, `Sz = -Re(Ey conj(Hx))/2`.

The exported signed and phase fields retain the complex solution. Display
normalization is recorded separately from the power accounting.

## Ordered implementation gate

1. **Completed audit and reliability foundation:** architecture, deployment,
   resource preflight, aggregate usage statistics, and responsive layouts.
2. **Completed:** shared finite BOX/handle geometry, legacy and cell-averaged
   scalar rasterization, and complex field/phase/relative Poynting export.
3. **Completed:** incoming Gaussian angular-spectrum equivalent current,
   homogeneous incident-power calibration, and incident/scattered separation.
4. **Completed with a stated phase limitation:** independent bidirectional
   solves and coefficient-magnitude/efficiency reciprocity. **Next:** common
   reference-plane phase de-embedding plus source-plane and window convergence.
5. **Then:** apodized/imported tooth tables, full-vector port contract and a
   separately reviewed full-vector finite propagator.

The remaining full-vector stage is not represented by inactive plots or copied
scalar results. It becomes user-facing only after its acceptance tests pass.
