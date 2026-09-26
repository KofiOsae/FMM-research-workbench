# Physics validation record

This workbench is a research preview. Numerical agreement on reference cases does not establish accuracy for every geometry. Record the exact material source, wavelength, Fourier basis, grid, and convergence study with each result.

## Independent checks

Run `.\.venv\Scripts\python.exe -m unittest -v` from this folder. The 55 tests currently pass.

| Calculation | Independent reference | Checked outcome |
| --- | --- | --- |
| Unpatterned interface | Fresnel s and p formulas at 0° and 30° | R and T agree to 10 decimal places; R+T=1 for lossless media. |
| Two-layer stack reversed | Reciprocity with incident/exit media exchanged | Total transmitted power agrees to 10 decimal places. |
| Uniform N-BK7 film | Independent Fabry–Pérot amplitude formula | Reflectance agrees to 10 decimal places; R+T=1. |
| Uniform absorbing gold film | Complex-index Fabry–Pérot formula | R and T agree to 10 decimal places at 0.905 µm; A is positive. |
| Vertical fields, material energy, and absorption | Fourier-domain unit-cell Poynting-flux loss versus global `1−R−T`; weak-loss derivative check; gold-film volume-loss integral | Complex E/H arrays are finite; summed layer absorption agrees with global A to 10⁻⁷; resolved local gold loss integrates to global absorption within 10⁻³; silica passes and gold fails the documented weak-loss energy gate. |
| Independent uniform-stack transfer matrix | FMM on a gold/silica two-film stack | R, T, and A agree within 10⁻⁸ for s, p, and unpolarized light at 0°, 35°, and 60°. Patterned layers are rejected by the transfer-matrix solver. |
| Oblique finite cell | Lossless disk on a 60° triangular/hexagonal Bravais cell | FMM conserves incident power to 10⁻⁷. |
| Rectangular and hexagonal bands | Reciprocal-lattice construction and symmetry paths | Finite TE/TM eigenvalues on Γ–X–M–Γ and Γ–M–K–Γ; monotonically accumulated path distance. |
| Measured-spectrum fit | Synthetic uniform-film spectrum with a hidden 0.26 µm thickness | Bounded fit recovers 0.26 µm within 10⁻⁵ and RMSE below 10⁻⁷. |
| Multi-parameter fit | Synthetic film with hidden thickness and refractive index, with reserved validation wavelengths | Both parameters recover within 10⁻⁵; training and validation RMSE are below 10⁻⁷; a 2×2 covariance/correlation report is produced. |
| Diffraction polarization channels | Isotropic interface under pure s and p illumination | Local p/s projections give S₁/S₀ = −1 for s and +1 for p, with S₂/S₀ = S₃/S₀ = 0, in both reflected and transmitted zeroth orders. |
| Far-field Jones/Stokes k-space map | Isotropic interface under coherent pure s and p illumination | Every valid transmitted sample retains S₁/S₀ = −1 for s or +1 for p; normalized Jones magnitude is one and unpolarized input is rejected. |
| Slab/grating phase matching | A slab mode with a grating period constructed as `2π/β` at normal incidence | The first reciprocal order gives zero relative momentum mismatch; uniform-layer effective index is recovered exactly. |
| Cavity Purcell estimate | Analytic single-mode Q/V expression | On-resonance unit-overlap result equals 3Q/(4π²) when V=(λ/n)³; detuning and overlap reduce it. |
| Planar dipole LDOS | Dyadic Green-function Sommerfeld integral | Index matching gives Γ/Γ₀=1 for both orientations and at every wavelength in a spectrum; a nearby lossy gold film gives a converged evanescent enhancement; isotropic averaging is (Γ⊥+2Γ∥)/3. |
| Planar objective collection | Angular-spectrum integration over the upper propagating cone | Collection is bounded by the upper-medium radiative power, respects `NA ≤ n`, and converges under quadrature refinement. It is collection efficiency, not a named-mode beta factor. |
| Full-vector rectangular waveguide | Fallahkhair Hx/Hy finite-difference eigenproblem | All six complex field components are finite; the silicon-core reference returns guided effective indices between cladding and core indices. Mesh, window, and boundary refinement remain required per design. |
| Bayesian spectrum inference | Synthetic bounded posterior with multiple independent chains | Posterior samples respect parameter bounds and report acceptance, split-chain R-hat, and effective sample size. The present likelihood is independent, constant-variance Gaussian. |
| Advanced local constitutive response | Homogeneous uniaxial Maxwell eigenproblem and analytic material conversions | Isotropic input returns degenerate bulk modes; rotations preserve tensor eigenvalues; thermal, Kerr, gain, and permeability conversions follow the displayed conventions. This is not a coupled stack solve. |
| Polarization winding arithmetic | Closed-loop doubled-angle unwrapping for a headless polarization axis | A synthetic radial vector texture returns charge +1 and a uniform texture returns zero. Structure-specific certification additionally requires every displayed physical and numerical check. |
| Bounded geometry optimization | Differential evolution with optional local polishing and linear constraints | A hidden uniform-film thickness is recovered from two wavelength targets within 0.015 µm and the selected design passes Fourier and grid validation. |
| Joint multi-resonance fitting and branch assignment | Synthetic spectrum containing one peak and one dip; reordered sweep features with supplied field overlap | BIC selects two jointly fitted components, both centers are recovered within 0.001 µm, and overlap-aware assignment preserves the two branches after their result order is reversed. |
| One-dimensional stripe | Grid is invariant along y | All y columns of permittivity agree; R+T=1 for a lossless case. |
| Unpolarized incidence | Separate s and p solves | Each reported power fraction is their arithmetic mean. |
| Uniform 2D crystal | Analytically folded plane waves | First TE/TM frequencies agree at sampled symmetry points. |
| Dielectric disk band | [MIT MPB tutorial](https://mpb.readthedocs.io/en/stable/Python_Tutorial/) | First TE frequency 0.37555 at the specified test point versus 0.372604 in MPB, a 0.79% difference at Fourier radius 7. |
| Symmetric planar waveguide | Independent even-mode half-slab equations `u tan(u)=w` for TE and `u tan(u)=(n_core/n_clad)² w` for TM | Effective indices agree to 11 decimal places in the tested case. |

The guided-mode solver additionally checks that the effective index lies strictly between the maximum cladding index and core index, and that its dispersion-equation residual is small. It reports how much an estimated group index changes when the wavelength differentiation step is halved. The core guided-power fraction uses |Eᵧ|² for TE and |Hᵧ|²/ε for TM.

## Before making a research claim

1. Check that every material is inside its cited wavelength range and appropriate for the specimen and temperature. The N-BK7, CaF₂, and fused-silica Sellmeier models return real n and omit absorption.
2. For spectra, inspect Fourier-order changes across three **actual** retained bases. At a wavelength you plan to quote, run the single-wavelength calculation and double the geometry grid. A passing 0.01 absolute-power-change screen is necessary for using this preview, but it is not an error bar.
3. For fields, inspect order and grid sensitivity of the map. A single bright pixel is not sufficient evidence of a physical hot spot.
4. For bands, compare each band across Fourier radii and check a denser Brillouin-zone path or independent solver before identifying a full gap.
5. For a waveguide, distinguish a bound in-plane eigenmode from an externally illuminated stack. The present model does not calculate coupling efficiency or guided-mode scattering at patterned interfaces.
6. Compare a case-specific result to an independent solver or measurement and include material and fabrication uncertainty before publication.

The historical patterned gold rectangle near 0.9 µm still fails its Fourier convergence check at the practical basis sizes in this app. Its spectra and fields should remain exploratory.






## Research workflow regression cases

| Workflow | Automated evidence | Scope |
|---|---|---|
| Complete-stack planar modes | TE and TM effective indices agree with the independent analytic symmetric-slab solver within 0.002 | Lossless laterally homogenized stacks |
| Fabrication tolerance | Seeded sampling, statistics, bounds, and feature tracking are deterministic | Independent Gaussian inputs |
| Band eigenmode fields | Uniform-medium reconstruction is finite, normalized, and has the expected qualified inversion overlap | Scalar 2D TE/TM PWE |
| Resonance field report | Four requested spectral positions are sampled and returned with R/T/A and field summaries | Driven vertical FMM fields |
| Experimental nuisance fit | Hidden wavelength offset, scale, and background are recovered from synthetic data | Conditional least-squares identification |
| Resonant Jones/Stokes subtraction | Coherent channel sidebands return finite residual Stokes data and a spacing-stability metric | Local background subtraction, not a pole residue |
| Uniform-stack complex pole | The outgoing transfer-matrix denominator converges below 1e-7 and a lossless slab gives matching total/radiative Q | Normal incidence, frozen indices, uniform layers |
| Methods report | Model, materials, versions, and reproducibility checklist are present | Documentation completeness |

| Vector-field arrays | Horizontal and vertical real E/H and Poynting arrays are finite and shape-consistent | Arrows are normalized display overlays; arrays remain quantitative |
| Expanded materials | Published wavelength nodes, range rejection, and passive-metal n+ik are regression tested | Isotropic bulk or cited-film models |

