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
TE finite-difference frequency-domain calculation with two independent source
choices on the same physical model. **Waveguide** launches the normalized
fundamental port mode and projects outgoing radiation on a Gaussian port.
**Gaussian** constructs the time-reversed angular spectrum, converts it to the
discrete equivalent current \(b=A_\mathrm{clad}E_\mathrm{inc}\), solves the
device, and projects the result on the right-going waveguide eigenmode. A
separate homogeneous solve supplies incident-power normalization and the
incident/scattered field separation.
The domain contains a finite grating, either a homogeneous lower half-space or
a finite BOX and handle, a uniform waveguide port, graded
complex-coordinate absorbers, and normalized modal and radiation monitors.
Its main target observable is

\[
\eta_{\mathrm{target}}=\frac{|\langle u_{\mathrm{target}},
E_{\mathrm{radiated}}\rangle_P|^2}{P_{\mathrm{target}}P_{\mathrm{incident}}}.
\]

The current one-mode scalar implementation reports guided output in each
direction, free-space reflection or upward radiation, substrate radiation, and
numerical residual. For passive lossless media the normalized budget should
satisfy

\[
1\approx P_\mathrm{guided,fwd} + P_\mathrm{guided,back}
+P_\mathrm{up}+P_\mathrm{sub}+\epsilon_\mathrm{num}.
\]

The interface reports target-mode efficiency, insertion loss, upward and
substrate radiation, residual forward guided power, back-reflection,
directionality, accounted power, and numerical/absorber residual. Its
validation command independently changes mesh spacing, absorber strength, and
all domain paddings, and interface rasterization. Its Trust table records the
declared tooth width, etch depth and removed area alongside the effective fill,
removed area, and nearest grid interfaces actually passed to the operator. This
separates a geometry-raster change from electromagnetic convergence. The selectable cell-averaged
epsilon corrects the scalar TE mass coefficient; it is not the anisotropic
subpixel tensor of a full-vector Maxwell discretization. Complex Ey, its phase,
relative reconstructed Poynting components, and the epsilon raster are exported.
A publication claim also requires wavelength sampling,
monitor-position checks, and comparison with an independent solver or
published benchmark.

This initial solver is restricted to reciprocal, isotropic, lossless 2D TE
structures invariant across their width. It does not yet support full-vector
3D focusing gratings, arbitrary finite-width masks, material absorption, or
several output modes. The reciprocity certificate now runs both source
constructions and compares their power-normalized complex coefficients. Its
current gate uses coefficient magnitude and efficiency. Raw phases are
exported, but phase reciprocity is not certified until both port reference
planes are de-embedded to a common origin. Those limitations are explicit in
every result export.

## 5. Bidirectional worked workflows

### Tutorial A — waveguide → grating → Gaussian port

1. Select **Waveguide mode → grating → Gaussian port**.
2. Define the SOI core, upper cladding, finite BOX and handle. The scalar port
   eigenproblem reports \(n_\mathrm{eff}\), propagation constant, core
   localization and normalization.
3. Estimate the first period from phase matching,
   \(\Lambda\approx\lambda_0/(n_\mathrm{eff}-n_c\sin\theta)\), then enter the
   actual finite period, fill factor, etch depth and number of teeth.
4. Define the outgoing Gaussian by its 1/e electric-field waist, center,
   signed global angle and reference phase.
5. Run the baseline controls. The analytic slab and nearly uniform guide are
   implementation controls; neither validates the patterned coupler.
6. Calculate the device. Inspect total complex field, phase and Poynting flow.
   Read Gaussian-port coupling together with upward radiation, substrate loss,
   guided reflection and the power residual.
7. Run spectrum, bounded optimization and fabrication tolerance. Treat the
   power residual, directionality and reflection as constraints.
8. Run Trust validation and the bidirectional reciprocity certificate. Export
   inputs, software settings, both complex coefficients and every certificate.

### Tutorial B — Gaussian port → grating → waveguide

1. Keep the same geometry and select **Gaussian port → grating → waveguide**.
2. The entered outgoing angle defines the port. The reciprocal incoming
   wavevector reverses sign in global coordinates; the interface reports both
   requested and actually sampled incoming angles.
3. Calculate the device. The homogeneous reference field normalizes incident
   power. Use the field selector to compare incident, scattered and total
   intensity or signed field.
4. Read **guided-mode coupling** as the right-going modal projection. Also
   inspect left-going guided power, reflected free-space power, substrate
   radiation and unaccounted/absorber residual.
5. Screen waist, position, angle and grating geometry. The optimizer and
   tolerance engine use the selected source direction and the same
   `selected_mode_coupling` observable.
6. Run Trust validation, then the reciprocity certificate. A failed
   reciprocity gate is evidence to refine mesh, window, absorbers, padding or
   reference-plane treatment; it must not be averaged away.

Both tutorials describe a 2D width-invariant Gaussian sheet and scalar TE
field. A realistic circular fiber above a finite-width focusing grating needs a
full-vector 3D solver and a two-dimensional transverse fiber-mode overlap.

## 6. Implemented finite-device design workflow

The finite-device workspace now carries one model through the complete nominal
design cycle:

1. **Baseline controls:** compare the finite-difference port effective index
   with the independent analytic asymmetric-slab dispersion relation and run
   the full domain in the nearly uniform-waveguide limit.
2. **Spectrum:** calculate target-mode efficiency versus wavelength and report
   sampled contiguous 1 dB and 3 dB bandwidths. Refine the wavelength grid
   around the selected peak before quoting bandwidth.
3. **Physical screens:** vary period, fill factor, etch depth, period count,
   target angle, or beam waist. A point is eligible only when its power-budget
   residual passes the stated screen.
4. **Bounded optimization:** vary one to five physical parameters together.
   Optimize the selected direction, the mean of both independently solved
   directions, or the weaker-direction efficiency. Minimum directionality,
   maximum reflection, and maximum power residual remain separate feasibility
   conditions and apply to both solves for a bidirectional objective.
5. **Robust optimization:** reuse fixed seeded Gaussian perturbations at every
   candidate and penalize variability. This is a screening objective, not a
   yield guarantee.
6. **Fabrication tolerance:** propagate stated independent Gaussian parameter
   distributions through the selected or both source directions and report
   objective percentiles plus yield against the entered efficiency,
   directionality, reflection, and residual limits.
7. **Trust validation:** independently vary mesh, absorber strength, all domain
   paddings, and binary versus cell-averaged scalar interfaces. Recompute a
   locally refined spectrum on the selected final design.

Project files and design-study exports preserve the finite-device inputs,
baseline controls, spectrum, parameter sweep, optimizer history, tolerance
samples, and Trust certificate.

## 7. Design variables and research outputs

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

## 8. Literature comparison case and scope gate

Marchetti et al., *Scientific Reports* **7**, 16670 (2017), provide a useful
high-efficiency evidence workflow: a full-vector 2D FDTD design, simultaneous
fill-factor and period apodization, etch optimization, a 260 nm SOI core,
oxide overcladding/BOX, a silicon handle, a 10.4 µm fiber mode-field diameter,
and fabrication measurements. They report 83% simulated peak coupling, 81%
best measured coupling, and a simulated 1 dB bandwidth of 32.8 nm.

Those headline values are **not** a numeric acceptance target for the current
scalar model. The **Load published 2017 tooth table** action imports all 24
published trench/tooth pairs, including the terminal trench, and sets the
reported 260 nm silicon layer, 160 nm etch, 2 µm BOX, 10° angle in top oxide,
and 5.4 µm 1/e field radius. The exact source data and scope notes are also in
`examples/marchetti_2017_apodized_scope_case.json`.

This is a geometry and workflow comparison. The paper used full-vector 2D FDTD,
a finite 680 nm top oxide followed by air, a 12 nm minimum conformal mesh, and a
finite 12 µm device width for the experimental fiber overlap. The Workbench
uses a scalar invariant-width operator and currently represents the top oxide
as a semi-infinite upper cladding. A direct numeric comparison with 83% would
therefore compare different electromagnetic models. Use the imported case to
verify geometry ingestion, raster error, power accounting, convergence sequence,
and trends; use a matched full-vector solver for quantitative reproduction.

## 9. Literature basis

- A. W. Snyder and J. D. Love, *Optical Waveguide Theory*, Chapman and Hall
  (1983): guided-mode orthogonality, excitation, and reciprocity.
- A. B. Fallahkhair, K. S. Li, and T. E. Murphy, “Vector finite difference
  modesolver for anisotropic dielectric waveguides,” *J. Lightwave Technol.*
  **26**, 1423–1431 (2008),
  [DOI 10.1109/JLT.2008.923643](https://doi.org/10.1109/JLT.2008.923643).
- D. Taillaert et al., “Grating couplers for coupling between optical fibers
  and nanophotonic waveguides,” *Jpn. J. Appl. Phys.* **45**, 6071–6077 (2006),
  [DOI 10.1143/JJAP.45.6071](https://doi.org/10.1143/JJAP.45.6071).
- R. Marchetti et al., “High-efficiency grating-couplers: demonstration of a
  new design strategy,” *Scientific Reports* **7**, 16670 (2017),
  [DOI 10.1038/s41598-017-16505-z](https://doi.org/10.1038/s41598-017-16505-z):
  simultaneous period/fill apodization, numerical optimization, tolerance,
  bandwidth, and experimental comparison.
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
