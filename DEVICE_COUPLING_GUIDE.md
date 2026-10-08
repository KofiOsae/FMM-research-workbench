# Research workflow for source, grating, and waveguide coupling

## 1. Distinguish the three electromagnetic quantities

The Workbench keeps three calculations separate because they answer different
questions and require different boundary conditions.

1. **Periodic diffraction efficiency** is obtained from the finite-stack
   Fourier modal method (FMM).  It reports power in plane-wave orders of an
   infinitely repeated unit cell.
2. **Uniform-port modal projection** expands a supplied complex field on the
   eigenmodes of a uniform waveguide cross section.  It reports source/mode
   compatibility at that plane.
3. **Finite-device coupling efficiency** must propagate a finite illumination
   field through a finite grating and project the resulting output field on
   normalized waveguide ports.

Neither (1) nor (2), nor their product, is generally equal to (3).  A product
is valid only after a common normalization and a derivation showing that all
multiple scattering, radiation, substrate, reflection, and transition channels
are represented.

## 2. Uniform-port formalism

Let the transverse port plane have normal \(\hat z\).  For a supplied field
\((\mathbf E_s,\mathbf H_s)\) and a normalized mode
\((\mathbf E_m,\mathbf H_m)\), the interface evaluates the electric-profile
diagnostic

\[
\eta_{E,m}=\frac{|\int_A \mathbf E_s\cdot\mathbf E_m^*\,dA|^2}
{\int_A|\mathbf E_s|^2dA\int_A|\mathbf E_m|^2dA}.
\]

When consistent E and H fields are supplied, the forward and backward
coefficients are evaluated by the Lorentz-reciprocity combinations of the
transverse fields.  The implementation uses

\[
a_m^+=\frac{I_1+I_2}{4\sqrt{P_sP_m}},\qquad
a_m^-=\frac{I_1-I_2}{4\sqrt{P_sP_m}},
\]

with

\[
I_1=\int_A( E_{s,x}H_{m,y}^*-E_{s,y}H_{m,x}^*)dA,
\]

\[
I_2=\int_A( E_{m,x}^*H_{s,y}-E_{m,y}^*H_{s,x})dA,
\]

and \(P=(1/2)\operatorname{Re}\int_A( E_xH_y^*-E_yH_x^*)dA\).
The reported directional efficiencies are \(|a_m^\pm|^2\).

### Validation gates

- Re-importing a solved mode on the same grid must return unit forward overlap
  and negligible backward overlap.
- Halving the mesh must stabilize both \(n_\mathrm{eff}\) and the selected
  overlap within stated tolerances.
- Domain padding must be varied independently of mesh spacing.
- The same physical mode must be tracked using effective index and field
  character; a list index alone is insufficient near crossings.
- Imported fields must cover a complete rectangular plane.  Values outside
  their measured/simulated rectangle are explicitly set to zero.
- The sum captured by requested bound modes is reported separately from the
  unresolved fraction.  The latter can contain radiation, continuum content,
  unrequested modes, and discretization error.

## 3. Current interface workflow

1. In **Guided modes**, define the waveguide, wavelength, material indices,
   domain padding, mesh, and requested modes.
2. Solve and identify the physical output mode using field shape,
   \(n_\mathrm{eff}\), TE-like fraction, and core localization.
3. Open **Uniform-port source → guided-mode projection**.
4. Use an analytic Gaussian/super-Gaussian or import complex field data.
5. Compare the normalized source, mode, and local overlap-density maps.
6. Screen waist, lateral position, tilt, or polarization.  The mode basis is
   solved once and reused because the waveguide and wavelength do not change.
7. Run **Validate mesh ×2** and then perform a separate padding comparison.
8. Save the project and export the projection JSON.  It includes the equations,
   mode table, solved-mode accounting, source screening, and mesh certificate.

## 4. Finite-device formulation implemented now

The **Finite grating coupler** workspace implements a two-dimensional scalar
TE finite-difference frequency-domain calculation. It launches the normalized
fundamental output-waveguide mode toward a finite, partially or fully etched
grating. Radiation into a selected free-space Gaussian mode is converted to
the reverse illumination-to-waveguide efficiency by Lorentz reciprocity. The
domain contains a finite grating, substrate, uniform output waveguide, graded
complex-coordinate absorbers, and normalized modal and radiation monitors.
Its main observable is

\[
\eta_{\mathrm{TE0}}^+=\frac{P_{\mathrm{TE0}}^+}{P_{\mathrm{source,in}}}.
\]

It must independently report \(\eta_{\mathrm{TE0}}^-\), reflected source
power, upward radiation, substrate radiation, material absorption, other bound
modes, and numerical residual.  For passive media the normalized budget should
satisfy

\[
1\approx \eta_{\mathrm{TE0}}^+ + \eta_{\mathrm{TE0}}^- + R
+P_\mathrm{up}+P_\mathrm{sub}+A+P_\mathrm{other}+\epsilon_\mathrm{num}.
\]

The interface reports target-mode efficiency, insertion loss, upward and
substrate radiation, residual forward guided power, back-reflection,
directionality, accounted power, and numerical/absorber residual. Its
validation command independently changes mesh spacing, absorber strength, and
all domain paddings. A publication claim also requires wavelength sampling,
monitor-position checks, and comparison with an independent solver or
published benchmark.

This initial solver is restricted to reciprocal, isotropic, lossless 2D TE
structures invariant across their width. It does not yet support full-vector
3D focusing gratings, arbitrary finite-width masks, material absorption, or
several output modes. Those limitations are explicit in every result export.

## 5. Design variables and research outputs

A finite-device design study should permit period, duty cycle, etch depth,
number of periods, chirp, apodization, fiber waist, fiber angle, lateral offset,
and output-waveguide geometry.  It should use the Workbench's existing reusable
observables and constraints, but rank candidates only after convergence-aware
reevaluation.  Nominal optimization should be followed by fabrication-aware
yield analysis using measured parameter distributions.

Every export should contain geometry and material provenance, source
normalization, port normalization, solver version, mesh/PML/domain convergence,
power budget, mode identity, objective and constraints, and a model
fingerprint.

## 6. Literature basis

- A. W. Snyder and J. D. Love, *Optical Waveguide Theory*, Chapman and Hall
  (1983): guided-mode orthogonality, excitation, and reciprocity.
- A. B. Fallahkhair, K. S. Li, and T. E. Murphy, “Vector finite difference
  modesolver for anisotropic dielectric waveguides,” *J. Lightwave Technol.*
  **26**, 1423–1431 (2008),
  [DOI 10.1109/JLT.2008.923643](https://doi.org/10.1109/JLT.2008.923643).
- D. Taillaert et al., “Grating couplers for coupling between optical fibers
  and nanophotonic waveguides,” *Jpn. J. Appl. Phys.* **45**, 6071–6077 (2006),
  [DOI 10.1143/JJAP.45.6071](https://doi.org/10.1143/JJAP.45.6071).
- A. F. Oskooi et al., “MEEP: A flexible free-software package for
  electromagnetic simulations by the FDTD method,” *Comput. Phys. Commun.*
  **181**, 687–702 (2010),
  [DOI 10.1016/j.cpc.2009.11.008](https://doi.org/10.1016/j.cpc.2009.11.008):
  finite-domain Maxwell solution and absorbing-boundary verification.
- A. Y. Piggott et al., “Inverse design and demonstration of a compact and
  broadband on-chip wavelength demultiplexer,” *Nature Photonics* **9**,
  374–377 (2015),
  [DOI 10.1038/nphoton.2015.69](https://doi.org/10.1038/nphoton.2015.69):
  constrained photonic-device optimization with fabrication considerations.

These sources establish the numerical and physical formalism.  They do not
validate a new device automatically; each claimed geometry requires its own
convergence, power-balance, reciprocity, and independent-comparison evidence.
