"""Two-dimensional finite grating out-coupling by scalar TE FDFD.

The implemented source launches a normalized bound waveguide mode toward a
finite grating.  A target-beam overlap is reported for that out-going field.
Lorentz reciprocity motivates a reverse-coupling comparison, but a genuine
free-space incident solve is deliberately not claimed by this module.

The time convention is exp(-i omega t).  With E = y-hat Ey, the solved equation
is (d_x^2 + d_z^2 + k0^2 epsilon_r) Ey = source for non-magnetic, isotropic
media.  Ey and its normal derivative are continuous across material interfaces.
This is a scalar, invariant-width approximation rather than a full-vector 3D
fiber-coupler calculation.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, spsolve
from scipy.optimize import differential_evolution
import run_jobs
from waveguide import WaveguideModel, solve_waveguide
import os


@dataclass(frozen=True)
class FiniteGratingModel:
    wavelength_um: float = 1.55
    core_n: float = 3.48
    substrate_n: float = 1.444
    cladding_n: float = 1.0
    box_thickness_um: float | None = None
    handle_n: float | None = None
    waveguide_height_um: float = .22
    etch_depth_um: float = .07
    period_um: float = .63
    fill_factor: float = .50
    periods: int = 8
    left_padding_um: float = 1.5
    right_padding_um: float = 2.0
    top_padding_um: float = 1.5
    substrate_depth_um: float = 1.5
    mesh_um: float = .05
    absorber_um: float = .50
    absorber_strength: float = 3.0
    discretization: str = "binary"
    subpixel_samples: int = 5
    target_angle_deg: float = -30.0
    target_waist_um: float = 3.0
    target_center_um: float | None = None
    target_phase_deg: float = 0.0

    def validate(self) -> None:
        values = [v for v in asdict(self).values() if isinstance(v, (int, float)) and v is not None]
        if not np.isfinite(values).all():
            raise ValueError("Finite-grating settings must be finite")
        if not .2 <= self.wavelength_um <= 20:
            raise ValueError("Wavelength must be 0.2–20 µm")
        if min(self.core_n, self.substrate_n, self.cladding_n) <= 0:
            raise ValueError("All refractive indices must be positive")
        if self.handle_n is not None and self.handle_n <= 0:
            raise ValueError("The silicon-handle refractive index must be positive")
        if not 1 <= self.periods <= 80 or not 0 < self.fill_factor < 1:
            raise ValueError("Use 1–80 periods and a fill factor strictly between zero and one")
        if not 0 < self.etch_depth_um <= self.waveguide_height_um:
            raise ValueError("Etch depth must be positive and no larger than the waveguide height")
        if not .005 <= self.mesh_um <= .15:
            raise ValueError("Finite-device mesh must be 0.005–0.15 µm; the cell-budget preflight may require a coarser value for a large domain")
        if min(self.left_padding_um, self.right_padding_um, self.top_padding_um,
               self.substrate_depth_um, self.absorber_um) <= 0:
            raise ValueError("Domain padding and absorber thickness must be positive")
        if self.absorber_um >= min(self.left_padding_um, self.right_padding_um,
                                   self.top_padding_um, self.substrate_depth_um):
            raise ValueError("Absorber thickness must be smaller than every domain padding")
        if self.box_thickness_um is not None:
            if self.box_thickness_um <= 0:
                raise ValueError("Finite BOX thickness must be positive")
            if self.handle_n is None:
                raise ValueError("A finite BOX requires the handle-substrate refractive index")
            clearance = self.box_thickness_um + self.absorber_um + 2*self.mesh_um
            if self.substrate_depth_um <= clearance:
                raise ValueError(
                    "Substrate depth must exceed BOX thickness + absorber thickness + two mesh cells so the lower monitor and absorber lie in the handle")
        if self.discretization not in ("binary", "cell_average"):
            raise ValueError("Discretization must be 'binary' or 'cell_average'")
        if not 2 <= self.subpixel_samples <= 9:
            raise ValueError("Use 2–9 subpixel samples per axis")
        if abs(self.target_angle_deg) >= 85 or self.target_waist_um <= 0:
            raise ValueError("Target angle must be below 85° and target waist must be positive")


def finite_grid_limit() -> int:
    """Return the configured hard cell limit for this service instance."""
    public = os.environ.get("PUBLIC_DEMO") == "1"
    default_limit = 80_000 if public else 300_000
    setting = "FINITE_GRID_LIMIT_PUBLIC" if public else "FINITE_GRID_LIMIT_LOCAL"
    try:
        return max(10_000, int(os.environ.get(setting, default_limit)))
    except ValueError:
        return default_limit


def _axes(model: FiniteGratingModel):
    length = model.periods*model.period_um
    x = np.arange(-model.left_padding_um, length+model.right_padding_um+model.mesh_um/2,
                  model.mesh_um)
    z = np.arange(-model.substrate_depth_um,
                  model.waveguide_height_um+model.top_padding_um+model.mesh_um/2,
                  model.mesh_um)
    public = os.environ.get("PUBLIC_DEMO") == "1"
    limit = finite_grid_limit()
    if x.size*z.size > limit:
        context = " on the shared online service" if public else ""
        minimum_mesh = np.sqrt((x[-1]-x[0])*(z[-1]-z[0])/limit)
        raise ValueError(
            f"Finite-device grid requests {x.size*z.size:,} cells{context}; the configured limit is {limit:,}. "
            f"For this domain use mesh about {minimum_mesh:.4f} µm or coarser, reduce padding/period count, or run the fine validation locally.")
    return x, z


def _point_permittivity(model: FiniteGratingModel, x, z, patterned=True):
    """Return exact-region epsilon at broadcast-compatible sample positions."""
    X, Z = np.broadcast_arrays(np.asarray(x), np.asarray(z))
    if model.box_thickness_um is None:
        below = np.full(X.shape, model.substrate_n**2)
    else:
        handle_n = model.handle_n if model.handle_n is not None else model.substrate_n
        below = np.where(Z >= -model.box_thickness_um,
                         model.substrate_n**2, handle_n**2)
    eps = np.where(Z < 0, below, model.cladding_n**2).astype(float)
    core = (Z >= 0) & (Z <= model.waveguide_height_um)
    eps[core] = model.core_n**2
    if patterned:
        length = model.periods*model.period_um
        phase = np.mod(X, model.period_um)/model.period_um
        grooves = ((X >= 0) & (X < length) &
                   (phase >= model.fill_factor))
        etched = ((Z > model.waveguide_height_um-model.etch_depth_um) &
                  (Z <= model.waveguide_height_um))
        eps[grooves & etched] = model.cladding_n**2
    return eps


def _permittivity(model: FiniteGratingModel, x, z):
    """Rasterize the device and its uniform port reference.

    ``binary`` reproduces the legacy node-sampled model. ``cell_average``
    averages epsilon over each finite-difference control cell by deterministic
    supersampling.  For this scalar TE equation epsilon is the mass coefficient,
    so volume averaging is the appropriate limited correction.  It is not the
    anisotropic subpixel tensor required by a full-vector discretization.
    """
    X, Z = np.meshgrid(x, z, indexing="ij")
    if model.discretization == "binary":
        eps = _point_permittivity(model, X, Z, True)
        reference_2d = _point_permittivity(model, X, Z, False)
    else:
        dx = float(x[1]-x[0]); dz = float(z[1]-z[0])
        offsets = (np.arange(model.subpixel_samples)+.5)/model.subpixel_samples-.5
        eps = np.zeros(X.shape, dtype=float)
        reference_2d = np.zeros(X.shape, dtype=float)
        for ox in offsets:
            for oz in offsets:
                eps += _point_permittivity(model, X+ox*dx, Z+oz*dz, True)
                reference_2d += _point_permittivity(model, X+ox*dx, Z+oz*dz, False)
        divisor = float(model.subpixel_samples**2)
        eps /= divisor; reference_2d /= divisor
    # The unpatterned reference is invariant along x apart from roundoff.
    return eps.astype(complex), reference_2d.mean(axis=0).astype(complex)


def _absorber_sigma(axis, low, high, thickness, strength):
    sigma = np.zeros_like(axis, dtype=float)
    lower = axis < low+thickness
    upper = axis > high-thickness
    sigma[lower] = strength*((low+thickness-axis[lower])/thickness)**3
    sigma[upper] = strength*((axis[upper]-(high-thickness))/thickness)**3
    return sigma


def _operator(eps, dx, dz, sx, sz, k0):
    """Assemble the five-point scalar Helmholtz operator directly as diagonals.

    Direct diagonal construction avoids the former per-entry Python row/column
    lists.  This materially lowers pre-factorization memory while preserving
    the same stencil and zero-Dirichlet outer rows.
    """
    nx, nz = eps.shape
    west=np.zeros((nx,nz),complex); east=np.zeros_like(west)
    south=np.zeros_like(west); north=np.zeros_like(west)
    interior=(slice(1,-1),slice(1,-1))
    sxm=(sx[1:-1]+sx[:-2])/2; sxp=(sx[1:-1]+sx[2:])/2
    szm=(sz[1:-1]+sz[:-2])/2; szp=(sz[1:-1]+sz[2:])/2
    west[interior]=1/(sx[1:-1,None]*sxm[:,None]*dx*dx)
    east[interior]=1/(sx[1:-1,None]*sxp[:,None]*dx*dx)
    south[interior]=1/(sz[None,1:-1]*szm[None,:]*dz*dz)
    north[interior]=1/(sz[None,1:-1]*szp[None,:]*dz*dz)
    main=k0*k0*eps-west-east-south-north
    boundary=np.zeros((nx,nz),bool); boundary[[0,-1],:]=True; boundary[:,[0,-1]]=True
    main[boundary]=1; west[boundary]=east[boundary]=south[boundary]=north[boundary]=0
    size=nx*nz
    return diags((west.ravel()[nz:],south.ravel()[1:],main.ravel(),
                  north.ravel()[:-1],east.ravel()[:-nz]),
                 (-nz,-1,0,1,nz),shape=(size,size),format="csr")


def _te_mode(reference_eps, z, k0):
    dz=float(z[1]-z[0]); interior=reference_eps[1:-1].real
    n=interior.size
    d2=diags([np.ones(n-1),-2*np.ones(n),np.ones(n-1)],[-1,0,1],format="csr")/dz**2
    matrix=d2+diags(k0*k0*interior,0)
    values,vectors=eigsh(matrix,k=min(16,max(1,n-2)),which="LA")
    order=np.argsort(values)[::-1]
    beta2=values[order]
    valid=np.where((beta2 > (k0*np.sqrt(min(reference_eps.real)))**2) &
                   (beta2 < (k0*np.sqrt(max(reference_eps.real))*1.001)**2))[0]
    if not valid.size:
        raise ValueError("No bound TE mode was found for the output waveguide")
    # A finite high-index handle also supports box-confined numerical modes.
    # Select the candidate localized in the device core above z=0 instead of
    # assuming that the largest propagation constant identifies the port.
    core_region=(z[1:-1]>=0) & (interior >= .99*max(reference_eps.real))
    scores=[]
    for candidate in valid:
        vector=vectors[:,order[candidate]]
        scores.append(float(np.sum(abs(vector[core_region])**2)/
                            max(np.sum(abs(vector)**2),1e-30)))
    selected=valid[int(np.argmax(scores))]
    if max(scores) < .05:
        raise ValueError("No device-layer-localized TE port mode was found; increase BOX thickness or the port-domain padding")
    beta=float(np.sqrt(beta2[selected])); phi=np.zeros(z.size,dtype=complex)
    phi[1:-1]=vectors[:,order[selected]]
    phi/=np.sqrt(np.trapezoid(abs(phi)**2,z))
    phase=phi[np.argmax(abs(phi))]
    phi*=np.exp(-1j*np.angle(phase))
    return beta,phi


def _mode_amplitudes(field, i, phi, z, dx, beta):
    e=field[i]; derivative=(field[i+1]-field[i-1])/(2*dx)
    c=np.trapezoid(np.conj(phi)*e,z)
    d=np.trapezoid(np.conj(phi)*derivative,z)
    return .5*(c+1j*d/beta), .5*(c-1j*d/beta)


def _angular_power(field, j, x, z, n, k0, upward=True):
    dz=float(z[1]-z[0]); dx=float(x[1]-x[0])
    derivative=(field[:,j+1]-field[:,j-1])/(2*dz)
    line=field[:,j]
    window=np.ones(x.size)
    edge=max(2,int(.08*x.size)); taper=.5*(1-np.cos(np.linspace(0,np.pi,edge)))
    window[:edge]=taper; window[-edge:]=taper[::-1]
    spectrum=dx*np.fft.fftshift(np.fft.fft(np.fft.ifftshift(line*window)))
    dspectrum=dx*np.fft.fftshift(np.fft.fft(np.fft.ifftshift(derivative*window)))
    kx=2*np.pi*np.fft.fftshift(np.fft.fftfreq(x.size,d=dx))
    kz2=(n*k0)**2-kx**2; propagating=kz2>0
    kz=np.sqrt(np.maximum(kz2,0))
    component=np.zeros_like(spectrum)
    if upward:
        component[propagating]=.5*(spectrum[propagating]+dspectrum[propagating]/(1j*kz[propagating]))
    else:
        component[propagating]=.5*(spectrum[propagating]-dspectrum[propagating]/(1j*kz[propagating]))
    dk=2*np.pi/(x.size*dx)
    power=float(np.sum(kz[propagating]*abs(component[propagating])**2)*dk/(2*np.pi))
    return kx,kz,component,propagating,power,dk


def solve_finite_grating(model: FiniteGratingModel) -> dict:
    model.validate(); x,z=_axes(model); dx=float(x[1]-x[0]); dz=float(z[1]-z[0])
    eps,reference=_permittivity(model,x,z); k0=2*np.pi/model.wavelength_um
    beta,phi=_te_mode(reference,z,k0)
    port_core_fraction=float(np.trapezoid(abs(phi[(z>=0)&(z<=model.waveguide_height_um)])**2,
                                          z[(z>=0)&(z<=model.waveguide_height_um)]))
    sigma_x=_absorber_sigma(x,x[0],x[-1],model.absorber_um,model.absorber_strength)
    sigma_z=_absorber_sigma(z,z[0],z[-1],model.absorber_um,model.absorber_strength)
    sx=1+1j*sigma_x; sz=1+1j*sigma_z
    actual_op=_operator(eps,dx,dz,sx,sz,k0)
    reference_eps=np.repeat(reference[None,:],x.size,axis=0)
    reference_op=_operator(reference_eps,dx,dz,sx,sz,k0)
    source_i=min(x.size-4, int(round((x[-1]-model.absorber_um-2*dx-x[0])/dx)))
    rhs=np.zeros_like(eps); rhs[source_i,:]=phi/dx
    rhs[[0,-1],:]=0; rhs[:,[0,-1]]=0
    reference_field=spsolve(reference_op,rhs.ravel()).reshape(eps.shape)
    total=spsolve(actual_op,rhs.ravel()).reshape(eps.shape)
    left_i=max(2,int(round((x[0]+model.absorber_um+3*dx-x[0])/dx)))
    right_i=max(left_i+3,source_i-4)
    leftward,rightward=_mode_amplitudes(total,left_i,phi,z,dx,beta)
    inc_left,inc_right=_mode_amplitudes(reference_field,right_i,phi,z,dx,beta)
    incident_power=beta*abs(inc_left)**2
    _,backward=_mode_amplitudes(total,right_i,phi,z,dx,beta)
    residual_guided=beta*abs(leftward)**2/incident_power
    back_reflection=beta*abs(backward)**2/incident_power
    top_j=min(z.size-3,int(round((z[-1]-model.absorber_um-3*dz-z[0])/dz)))
    bottom_j=max(2,int(round((z[0]+model.absorber_um+3*dz-z[0])/dz)))
    kx,kz,up,prop,up_power,dk=_angular_power(total,top_j,x,z,model.cladding_n,k0,True)
    lower_n = (model.handle_n if model.box_thickness_um is not None
               and model.handle_n is not None else model.substrate_n)
    *_,down_power,_=_angular_power(total,bottom_j,x,z,lower_n,k0,False)
    upward=up_power/incident_power; substrate=down_power/incident_power
    center=model.target_center_um if model.target_center_um is not None else model.periods*model.period_um/2
    target=np.exp(-((x-center)/model.target_waist_um)**2)*np.exp(
        1j*(model.cladding_n*k0*np.sin(np.deg2rad(model.target_angle_deg))*x+
            np.deg2rad(model.target_phase_deg)))
    target_spectrum=dx*np.fft.fftshift(np.fft.fft(np.fft.ifftshift(target)))
    target_norm=float(np.sum(kz[prop]*abs(target_spectrum[prop])**2)*dk/(2*np.pi))
    overlap=np.sum(kz[prop]*np.conj(target_spectrum[prop])*up[prop])*dk/(2*np.pi)
    target_eff=float(abs(overlap)**2/(target_norm*incident_power)) if target_norm else 0
    accounted=residual_guided+back_reflection+upward+substrate
    stride=max(1,int(max(x.size,z.size)/240))
    sampled=total[::stride,::stride]
    intensity=abs(sampled)**2; field_scale=max(float(np.sqrt(intensity.max())),1e-30)
    intensity/=field_scale**2
    d_ey_dx=np.gradient(total,dx,axis=0); d_ey_dz=np.gradient(total,dz,axis=1)
    hx=1j*d_ey_dz/k0; hz=-1j*d_ey_dx/k0
    sx_power=.5*np.real(total*np.conj(hz)); sz_power=-.5*np.real(total*np.conj(hx))
    power_scale=max(float(np.max(np.sqrt(sx_power**2+sz_power**2))),1e-30)
    sampled_eps=eps[::stride,::stride].real
    stack = ({"type":"semi_infinite_substrate","substrate_n":model.substrate_n}
             if model.box_thickness_um is None else
             {"type":"finite_box_and_handle","box_n":model.substrate_n,
              "box_thickness_um":model.box_thickness_um,"handle_n":model.handle_n})
    return {"method":"2D scalar TE finite-difference frequency domain; guided-mode out-coupling solve",
        "model":asdict(model),"grid":{"nx":x.size,"nz":z.size,"cells":x.size*z.size,"dx_um":dx,"dz_um":dz},
        "geometry":{"upper_cladding_n":model.cladding_n,"device_n":model.core_n,
                    "lower_stack":stack,"rasterization":model.discretization,
                    "subpixel_samples_per_axis":model.subpixel_samples if model.discretization=="cell_average" else 1},
        "formulation":{"time_convention":"exp(-i omega t)",
            "equation":"(d_x^2 + d_z^2 + k0^2 epsilon_r) Ey = b",
            "interface_conditions":"Ey and its normal derivative are continuous for this scalar TE, non-magnetic model",
            "magnetic_reconstruction":"Hx = i d_z(Ey)/k0; Hz = -i d_x(Ey)/k0 in relative units",
            "source":"discrete line source calibrated by an unpatterned reference-waveguide solve",
            "excitation":"waveguide_to_grating_to_target_mode"},
        "mode":{"n_eff":beta/k0,"beta_per_um":beta,
                "device_core_fraction":port_core_fraction,
                "selection":"maximum electric-profile localization in the device layer above z=0",
                "normalization":"integral |Ey|^2 dz = 1"},
        "efficiencies":{"target_free_space_mode":target_eff,"upward_radiation":upward,
            "substrate_radiation":substrate,"residual_forward_waveguide":residual_guided,
            "back_reflection":back_reflection,"accounted_power":accounted,
            "numerical_or_absorber_residual":1-accounted,
            "insertion_loss_db":float(-10*np.log10(max(target_eff,1e-30))),
            "directionality":float(upward/max(upward+substrate,1e-30))},
        "field":{"x_um":x[::stride].tolist(),"z_um":z[::stride].tolist(),
                 "normalization":"Ey is divided by max|Ey|; relative Poynting components are divided by max sqrt(Sx^2+Sz^2)",
                 "normalized_Ey2":intensity.T.tolist(),
                 "normalized_Ey_real":(sampled.real/field_scale).T.tolist(),
                 "normalized_Ey_imag":(sampled.imag/field_scale).T.tolist(),
                 "Ey_phase_rad":np.angle(sampled).T.tolist(),
                 "normalized_Sx":(sx_power[::stride,::stride]/power_scale).T.tolist(),
                 "normalized_Sz":(sz_power[::stride,::stride]/power_scale).T.tolist(),
                 "epsilon_r":sampled_eps.T.tolist()},
        "reciprocity_statement":"This result is one genuine waveguide-source solve. For reciprocal isotropic materials it defines the corresponding reciprocal-port overlap, but it is not an independently solved fiber-incident field and is not used as a reciprocity-error test.",
        "scope":"2D scalar TE finite-device out-coupling solver, invariant across width. It includes finite length, partial etch, optional finite BOX and handle, substrate leakage, back-reflection, absorbing boundaries, and target-mode overlap. It does not include finite lateral width, full-vector 3D polarization mixing, or a validated incoming-fiber source.",
        "references":[{"title":"Optical Waveguide Theory","authors":"Snyder and Love","applies_to":"mode normalization and reciprocity"},
          {"title":"Grating couplers for coupling between optical fibers and nanophotonic waveguides","doi":"10.1143/JJAP.45.6071"},
          {"title":"Improving accuracy by subpixel smoothing in the finite-difference time domain","doi":"10.1364/OL.31.002972","applies_to":"interface-aware discretization context; this scalar implementation uses cell-averaged epsilon, not the paper's full tensor method"},
          {"title":"MEEP: A flexible free-software package for electromagnetic simulations by the FDTD method","doi":"10.1016/j.cpc.2009.11.008","applies_to":"finite-domain validation context"}]}


def validate_finite_grating(model: FiniteGratingModel) -> dict:
    coarse=solve_finite_grating(model)
    refined=solve_finite_grating(FiniteGratingModel(**{**asdict(model),"mesh_um":model.mesh_um/1.5}))
    keys=("target_free_space_mode","upward_radiation","substrate_radiation","back_reflection")
    changes={key:abs(refined["efficiencies"][key]-coarse["efficiencies"][key]) for key in keys}
    absorber=solve_finite_grating(FiniteGratingModel(**{**asdict(model),
        "absorber_strength":model.absorber_strength*1.25}))
    absorber_changes={key:abs(absorber["efficiencies"][key]-coarse["efficiencies"][key]) for key in keys}
    padding=solve_finite_grating(FiniteGratingModel(**{**asdict(model),
        "left_padding_um":model.left_padding_um*1.15,
        "right_padding_um":model.right_padding_um*1.15,
        "top_padding_um":model.top_padding_um*1.15,
        "substrate_depth_um":model.substrate_depth_um*1.15}))
    padding_changes={key:abs(padding["efficiencies"][key]-coarse["efficiencies"][key]) for key in keys}
    alternate_discretization = "binary" if model.discretization == "cell_average" else "cell_average"
    interface=solve_finite_grating(FiniteGratingModel(**{**asdict(model),
        "discretization":alternate_discretization}))
    interface_changes={key:abs(interface["efficiencies"][key]-coarse["efficiencies"][key]) for key in keys}
    power_residuals={name:abs(run["efficiencies"]["numerical_or_absorber_residual"])
                     for name,run in (("coarse",coarse),("refined",refined),
                                      ("absorber",absorber),("padding",padding),
                                      ("interface",interface))}
    maximum=max([*changes.values(),*absorber_changes.values(),*padding_changes.values(),
                 *interface_changes.values()])
    research=maximum <= 1e-2 and max(power_residuals.values()) <= 1e-2
    publication=maximum <= 2e-3 and max(power_residuals.values()) <= 2e-3
    return {"coarse":coarse,"refined":refined,"absorber_variant":absorber,
            "padding_variant":padding,"interface_variant":interface,
            "interface_variant_name":alternate_discretization,"changes":changes,
            "absorber_changes":absorber_changes,"padding_changes":padding_changes,
            "interface_changes":interface_changes,
            "power_budget_residuals":power_residuals,"maximum_change":maximum,
            "passed_research":research,"passed_publication_screen":publication,
            "interpretation":"The certificate varies mesh, absorber strength, all domain paddings, and binary versus cell-averaged interface rasterization independently. Cell averaging is the scalar TE mass-term correction, not full-vector tensor smoothing. Wavelength sampling and comparison with an external benchmark remain separate publication checks."}


def benchmark_finite_grating(model: FiniteGratingModel) -> dict:
    """Run independent modal and uniform-device baseline controls.

    The analytic asymmetric-slab dispersion relation is independent of the
    finite-difference port eigensolver.  A nearly unetched uniform guide then
    exposes spurious radiation, reflection, or absorber error in the full
    finite-domain calculation.
    """
    analytic = solve_waveguide(WaveguideModel(
        wavelength_um=model.wavelength_um,
        thickness_um=model.waveguide_height_um,
        core_material="dielectric", core_n=model.core_n,
        top_material="dielectric", top_n=model.cladding_n,
        bottom_material="dielectric", bottom_n=model.substrate_n,
        polarization="TE"))
    if not analytic["modes"]:
        raise ValueError("The analytic slab control found no bound TE mode")
    uniform = solve_finite_grating(FiniteGratingModel(**{
        **asdict(model), "fill_factor": .999999,
        "target_center_um": None}))
    numerical_neff = float(uniform["mode"]["n_eff"])
    analytic_neff = float(analytic["modes"][0]["n_eff"])
    neff_error = abs(numerical_neff-analytic_neff)
    e = uniform["efficiencies"]
    spurious = e["upward_radiation"]+e["substrate_radiation"]+e["back_reflection"]
    # This is a deliberately loose baseline gate because the scalar port uses
    # the same Cartesian staircasing as the device mesh.  The separate Trust
    # certificate is the quantitative mesh-convergence gate.
    passed = neff_error <= .15 and spurious <= 2e-2 and abs(e["numerical_or_absorber_residual"]) <= 2e-2
    return {
        "analytic_slab_n_eff": analytic_neff,
        "finite_difference_port_n_eff": numerical_neff,
        "absolute_n_eff_difference": neff_error,
        "uniform_waveguide": _compact(uniform),
        "spurious_scattered_fraction": float(spurious),
        "passed": bool(passed),
        "criteria": {"maximum_absolute_n_eff_difference": .15,
                     "maximum_spurious_scattered_fraction": .02,
                     "maximum_power_budget_residual": .02},
        "interpretation": "This baseline checks the port mode against an independent analytic asymmetric-slab dispersion relation and checks the full finite domain in the uniform-waveguide limit. The n_eff gate is a coarse implementation check because the finite-difference port inherits Cartesian staircasing; use the Trust mesh study for quantitative convergence. This is not external validation of a patterned coupler.",
        "external_validation_required": "External validation remains required for a publication claim: compare the final patterned-device observable with a matched literature case, an independent full-wave solver, or experiment."
    }


_DESIGN_PARAMETERS = {name for name in FiniteGratingModel.__dataclass_fields__
                      if name not in {"target_center_um"}}
def _with(model: FiniteGratingModel, parameter: str, value: float) -> FiniteGratingModel:
    if parameter not in _DESIGN_PARAMETERS:
        raise ValueError(f"Unsupported finite-device parameter {parameter}")
    if parameter == "periods":
        value = int(round(value))
    return FiniteGratingModel(**{**asdict(model), parameter: value})


def _compact(result: dict) -> dict:
    return {"wavelength_um": result["model"]["wavelength_um"],
            "n_eff": result["mode"]["n_eff"], **result["efficiencies"]}


def finite_grating_spectrum(model: FiniteGratingModel, start: float, stop: float,
                            points: int) -> dict:
    if not np.isfinite([start,stop]).all() or start >= stop or not 3 <= points <= 81:
        raise ValueError("Use 3–81 wavelengths with increasing finite limits")
    rows=[]
    for index,wavelength in enumerate(np.linspace(start,stop,points)):
        run_jobs.progress(index,points,"Finite grating spectrum")
        rows.append(_compact(solve_finite_grating(_with(model,"wavelength_um",float(wavelength)))))
    values=np.array([row["target_free_space_mode"] for row in rows]); peak=int(np.argmax(values))
    def width(drop_db):
        threshold=values[peak]*10**(-drop_db/10); valid=values>=threshold
        lo=peak; hi=peak
        while lo>0 and valid[lo-1]: lo-=1
        while hi+1<len(rows) and valid[hi+1]: hi+=1
        return float(rows[hi]["wavelength_um"]-rows[lo]["wavelength_um"]), [lo,hi]
    one_db,one_bounds=width(1); three_db,three_bounds=width(3)
    return {"rows":rows,"peak":rows[peak],"peak_index":peak,
            "bandwidth_1db_um":one_db,"bandwidth_3db_um":three_db,
            "bandwidth_1db_indices":one_bounds,"bandwidth_3db_indices":three_bounds,
            "interpretation":"Bandwidth uses contiguous sampled points around the maximum. Refine wavelength sampling before publication."}


def finite_grating_sweep(model: FiniteGratingModel, parameter: str,
                         start: float, stop: float, points: int) -> dict:
    if parameter not in _DESIGN_PARAMETERS or parameter in {"mesh_um","absorber_um","absorber_strength"}:
        raise ValueError("Choose a physical grating, waveguide, or target-beam parameter")
    if not np.isfinite([start,stop]).all() or start >= stop or not 3 <= points <= 61:
        raise ValueError("Use 3–61 points with increasing finite limits")
    rows=[]
    for index,value in enumerate(np.linspace(start,stop,points)):
        run_jobs.progress(index,points,"Finite grating parameter sweep")
        result=solve_finite_grating(_with(model,parameter,float(value)))
        rows.append({"parameter_value":int(round(value)) if parameter=="periods" else float(value),
                     **_compact(result)})
    eligible=[row for row in rows if abs(row["numerical_or_absorber_residual"]) <= .03]
    if not eligible:
        raise ValueError("No sweep point met the 3% power-budget screen; refine the domain before optimization")
    best=max(eligible,key=lambda row:row["target_free_space_mode"])
    return {"parameter":parameter,"rows":rows,"best":best,
            "eligibility":"|numerical or absorber residual| <= 0.03"}


def finite_grating_tolerance(model: FiniteGratingModel, uncertainties: list[dict],
                             samples: int = 20, seed: int = 12345,
                             minimum_efficiency: float = .5,
                             maximum_reflection: float = .05,
                             minimum_directionality: float = .5) -> dict:
    if not 5 <= samples <= 200:
        raise ValueError("Use 5–200 tolerance samples")
    normalized=[]
    for item in uncertainties:
        name=str(item["parameter"]); sigma=float(item["sigma"])
        if name not in _DESIGN_PARAMETERS or name in {"mesh_um","absorber_um","absorber_strength"} or sigma <= 0:
            raise ValueError("Tolerance parameters must be physical finite-device inputs with positive sigma")
        lower=float(item.get("lower",-np.inf)); upper=float(item.get("upper",np.inf))
        if lower >= upper:
            raise ValueError("Every fabrication bound must have lower < upper")
        normalized.append({"parameter":name,"sigma":sigma,
            "lower":lower,"upper":upper})
    if not normalized or len(normalized)>6:
        raise ValueError("Provide 1–6 uncertain parameters")
    rng=np.random.default_rng(seed); rows=[]
    for sample in range(samples):
        run_jobs.progress(sample,samples,"Finite grating fabrication tolerance")
        values=asdict(model)
        for item in normalized:
            draw=float(np.clip(rng.normal(float(values[item["parameter"]]),item["sigma"]),item["lower"],item["upper"]))
            values[item["parameter"]]=int(round(draw)) if item["parameter"]=="periods" else draw
        try:
            result=solve_finite_grating(FiniteGratingModel(**values)); row={"sample":sample,"parameters":{i["parameter"]:values[i["parameter"]] for i in normalized},**_compact(result)}
            row["passed_specification"]=bool(row["target_free_space_mode"]>=minimum_efficiency and row["back_reflection"]<=maximum_reflection and row["directionality"]>=minimum_directionality and abs(row["numerical_or_absorber_residual"])<=.03)
            rows.append(row)
        except (ValueError,np.linalg.LinAlgError) as exc:
            rows.append({"sample":sample,"parameters":{i["parameter"]:values[i["parameter"]] for i in normalized},"error":str(exc),"passed_specification":False})
    valid=[row for row in rows if "error" not in row]
    if not valid: raise ValueError("Every fabrication sample failed")
    values=np.array([row["target_free_space_mode"] for row in valid])
    serialized_uncertainties=[{**item,
        "lower":item["lower"] if np.isfinite(item["lower"]) else None,
        "upper":item["upper"] if np.isfinite(item["upper"]) else None}
        for item in normalized]
    return {"samples":rows,"completed":len(valid),"failed":samples-len(valid),"seed":seed,
        "uncertainties":serialized_uncertainties,"specification":{"minimum_efficiency":minimum_efficiency,
        "maximum_reflection":maximum_reflection,"minimum_directionality":minimum_directionality,
        "maximum_power_residual":.03},"yield_fraction":float(sum(r["passed_specification"] for r in rows)/samples),
        "efficiency_statistics":{"mean":float(values.mean()),"sd":float(values.std(ddof=1)) if len(values)>1 else 0,
        "p05":float(np.percentile(values,5)),"median":float(np.median(values)),"p95":float(np.percentile(values,95))},
        "interpretation":"Yield is conditional on the entered independent Gaussian fabrication model and numerical power-budget screen."}


def optimize_finite_grating(model: FiniteGratingModel, variables: list[dict],
                            generations: int = 3, population: int = 5,
                            seed: int = 12345, minimum_directionality: float = 0,
                            maximum_reflection: float = 1,
                            maximum_power_residual: float = .03,
                            robust_samples: int = 1,
                            uncertainty_sigma: dict | None = None,
                            variability_weight: float = 0) -> dict:
    if not 1 <= len(variables) <= 5 or not 1 <= generations <= 30 or not 4 <= population <= 15:
        raise ValueError("Use 1–5 variables, 1–30 generations, and population 4–15")
    names=[]; bounds=[]
    for item in variables:
        name=str(item["parameter"]); lower=float(item["lower"]); upper=float(item["upper"])
        if name not in _DESIGN_PARAMETERS or name in {"mesh_um","absorber_um","absorber_strength"} or not lower<upper or name in names:
            raise ValueError("Optimization variables must be unique physical inputs with increasing bounds")
        names.append(name); bounds.append((lower,upper))
    if not 1 <= robust_samples <= 9: raise ValueError("Use 1–9 fixed robust samples per candidate")
    sigmas={str(k):float(v) for k,v in (uncertainty_sigma or {}).items() if str(k) in names and float(v)>0}
    rng=np.random.default_rng(seed); perturb=rng.normal(size=(robust_samples,len(names))); history=[]
    def evaluate(vector):
        losses=[]
        for sample in range(robust_samples):
            values=[]
            for i,(name,(lower,upper)) in enumerate(zip(names,bounds)):
                value=float(vector[i]+perturb[sample,i]*sigmas.get(name,0)); values.append(float(np.clip(value,lower,upper)))
            candidate=model
            for name,value in zip(names,values): candidate=_with(candidate,name,value)
            try:
                result=solve_finite_grating(candidate); e=result["efficiencies"]
                violation=max(0,minimum_directionality-e["directionality"])+max(0,e["back_reflection"]-maximum_reflection)+max(0,abs(e["numerical_or_absorber_residual"])-maximum_power_residual)
                losses.append(-e["target_free_space_mode"]+1e3*violation)
            except (ValueError,np.linalg.LinAlgError): losses.append(1e6)
        score=float(np.mean(losses)+variability_weight*np.std(losses)); history.append({"parameters":{k:float(v) for k,v in zip(names,vector)},"loss":score})
        run_jobs.progress(len(history),max(1,population*len(names)*(generations+1)),"Finite grating optimization")
        return score
    result=differential_evolution(evaluate,bounds,seed=seed,maxiter=generations,popsize=population,polish=False,updating="immediate",workers=1)
    best=model
    for name,value in zip(names,result.x): best=_with(best,name,float(value))
    solved=solve_finite_grating(best); e=solved["efficiencies"]
    nominal_feasible=e["directionality"]>=minimum_directionality and e["back_reflection"]<=maximum_reflection and abs(e["numerical_or_absorber_residual"])<=maximum_power_residual
    robust_efficiencies=[]; robust_pass=[]
    for sample in range(robust_samples):
        candidate=best
        for i,(name,(lower,upper)) in enumerate(zip(names,bounds)):
            value=float(np.clip(result.x[i]+perturb[sample,i]*sigmas.get(name,0),lower,upper))
            candidate=_with(candidate,name,value)
        trial=solve_finite_grating(candidate)["efficiencies"]
        robust_efficiencies.append(trial["target_free_space_mode"])
        robust_pass.append(trial["directionality"]>=minimum_directionality and
            trial["back_reflection"]<=maximum_reflection and
            abs(trial["numerical_or_absorber_residual"])<=maximum_power_residual)
    robust_summary={"mean_efficiency":float(np.mean(robust_efficiencies)),
        "sd_efficiency":float(np.std(robust_efficiencies)),
        "minimum_efficiency":float(np.min(robust_efficiencies)),
        "all_samples_feasible":bool(all(robust_pass))}
    feasible=nominal_feasible and robust_summary["all_samples_feasible"]
    return {"best_parameters":{k:(int(round(v)) if k=="periods" else float(v)) for k,v in zip(names,result.x)},
        "best_result":solved,"best_efficiency":e["target_free_space_mode"],"feasible":feasible,
        "constraints":{"minimum_directionality":minimum_directionality,"maximum_reflection":maximum_reflection,"maximum_power_residual":maximum_power_residual},
        "robust_samples":robust_samples,"uncertainty_sigma":sigmas,"variability_weight":variability_weight,
        "robust_summary":robust_summary,
        "evaluations":len(history),"history":history,"seed":seed,
        "warning":"Bounded stochastic optimization does not prove a global optimum. Run independent seeds, Trust validation, wavelength sampling, and fabrication tolerance before accepting the design."}
