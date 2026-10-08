"""Two-dimensional finite grating-to-waveguide coupling by scalar TE FDFD.

The reciprocal calculation launches a normalized bound waveguide mode toward a
finite grating.  Radiation into a requested free-space beam is equal to the
reverse beam-to-waveguide efficiency for reciprocal isotropic materials.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh, spsolve


@dataclass(frozen=True)
class FiniteGratingModel:
    wavelength_um: float = 1.55
    core_n: float = 3.48
    substrate_n: float = 1.444
    cladding_n: float = 1.0
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
    target_angle_deg: float = -30.0
    target_waist_um: float = 3.0
    target_center_um: float | None = None

    def validate(self) -> None:
        values = [v for v in asdict(self).values() if isinstance(v, (int, float)) and v is not None]
        if not np.isfinite(values).all():
            raise ValueError("Finite-grating settings must be finite")
        if not .2 <= self.wavelength_um <= 20:
            raise ValueError("Wavelength must be 0.2–20 µm")
        if not 1 <= self.periods <= 80 or not 0 < self.fill_factor < 1:
            raise ValueError("Use 1–80 periods and a fill factor strictly between zero and one")
        if not 0 < self.etch_depth_um <= self.waveguide_height_um:
            raise ValueError("Etch depth must be positive and no larger than the waveguide height")
        if not .01 <= self.mesh_um <= .15:
            raise ValueError("Finite-device mesh must be 0.01–0.15 µm")
        if min(self.left_padding_um, self.right_padding_um, self.top_padding_um,
               self.substrate_depth_um, self.absorber_um) <= 0:
            raise ValueError("Domain padding and absorber thickness must be positive")
        if self.absorber_um >= min(self.left_padding_um, self.right_padding_um,
                                   self.top_padding_um, self.substrate_depth_um):
            raise ValueError("Absorber thickness must be smaller than every domain padding")
        if abs(self.target_angle_deg) >= 85 or self.target_waist_um <= 0:
            raise ValueError("Target angle must be below 85° and target waist must be positive")


def _axes(model: FiniteGratingModel):
    length = model.periods*model.period_um
    x = np.arange(-model.left_padding_um, length+model.right_padding_um+model.mesh_um/2,
                  model.mesh_um)
    z = np.arange(-model.substrate_depth_um,
                  model.waveguide_height_um+model.top_padding_um+model.mesh_um/2,
                  model.mesh_um)
    if x.size*z.size > 120_000:
        raise ValueError(f"Finite-device grid requests {x.size*z.size:,} cells; reduce the domain or use a coarser mesh (limit 120,000)")
    return x, z


def _permittivity(model: FiniteGratingModel, x, z):
    eps = np.where(z[None, :] < 0, model.substrate_n**2, model.cladding_n**2)
    eps = np.repeat(eps, x.size, axis=0).astype(complex)
    core = (z >= 0) & (z <= model.waveguide_height_um)
    eps[:, core] = model.core_n**2
    length = model.periods*model.period_um
    phase = np.mod(x, model.period_um)/model.period_um
    grooves = (x >= 0) & (x < length) & (phase >= model.fill_factor)
    etched = (z > model.waveguide_height_um-model.etch_depth_um) & (z <= model.waveguide_height_um)
    eps[np.ix_(grooves, etched)] = model.cladding_n**2
    reference = np.where(z < 0, model.substrate_n**2, model.cladding_n**2).astype(complex)
    reference[core] = model.core_n**2
    return eps, reference


def _absorber_sigma(axis, low, high, thickness, strength):
    sigma = np.zeros_like(axis, dtype=float)
    lower = axis < low+thickness
    upper = axis > high-thickness
    sigma[lower] = strength*((low+thickness-axis[lower])/thickness)**3
    sigma[upper] = strength*((axis[upper]-(high-thickness))/thickness)**3
    return sigma


def _operator(eps, dx, dz, sx, sz, k0):
    nx, nz = eps.shape
    rows, cols, vals = [], [], []
    def add(r, c, v): rows.append(r); cols.append(c); vals.append(v)
    for i in range(nx):
        for j in range(nz):
            q = i*nz+j
            if i in (0, nx-1) or j in (0, nz-1):
                add(q, q, 1.0)
                continue
            sxm=(sx[i]+sx[i-1])/2; sxp=(sx[i]+sx[i+1])/2
            szm=(sz[j]+sz[j-1])/2; szp=(sz[j]+sz[j+1])/2
            axm=1/(sx[i]*sxm*dx*dx); axp=1/(sx[i]*sxp*dx*dx)
            azm=1/(sz[j]*szm*dz*dz); azp=1/(sz[j]*szp*dz*dz)
            add(q,(i-1)*nz+j,axm); add(q,(i+1)*nz+j,axp)
            add(q,i*nz+j-1,azm); add(q,i*nz+j+1,azp)
            add(q,q,k0*k0*eps[i,j]-axm-axp-azm-azp)
    return coo_matrix((vals,(rows,cols)),shape=(nx*nz,nx*nz)).tocsr()


def _te_mode(reference_eps, z, k0):
    dz=float(z[1]-z[0]); interior=reference_eps[1:-1].real
    n=interior.size
    d2=diags([np.ones(n-1),-2*np.ones(n),np.ones(n-1)],[-1,0,1],format="csr")/dz**2
    matrix=d2+diags(k0*k0*interior,0)
    values,vectors=eigsh(matrix,k=min(6,max(1,n-2)),which="LA")
    order=np.argsort(values)[::-1]
    beta2=values[order]
    valid=np.where((beta2 > (k0*np.sqrt(min(reference_eps.real)))**2) &
                   (beta2 < (k0*np.sqrt(max(reference_eps.real))*1.001)**2))[0]
    if not valid.size:
        raise ValueError("No bound TE mode was found for the output waveguide")
    beta=float(np.sqrt(beta2[valid[0]])); phi=np.zeros(z.size,dtype=complex)
    phi[1:-1]=vectors[:,order[valid[0]]]
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
    *_,down_power,_=_angular_power(total,bottom_j,x,z,model.substrate_n,k0,False)
    upward=up_power/incident_power; substrate=down_power/incident_power
    center=model.target_center_um if model.target_center_um is not None else model.periods*model.period_um/2
    target=np.exp(-((x-center)/model.target_waist_um)**2)*np.exp(1j*model.cladding_n*k0*np.sin(np.deg2rad(model.target_angle_deg))*x)
    target_spectrum=dx*np.fft.fftshift(np.fft.fft(np.fft.ifftshift(target)))
    target_norm=float(np.sum(kz[prop]*abs(target_spectrum[prop])**2)*dk/(2*np.pi))
    overlap=np.sum(kz[prop]*np.conj(target_spectrum[prop])*up[prop])*dk/(2*np.pi)
    target_eff=float(abs(overlap)**2/(target_norm*incident_power)) if target_norm else 0
    accounted=residual_guided+back_reflection+upward+substrate
    stride=max(1,int(max(x.size,z.size)/240))
    intensity=abs(total[::stride,::stride])**2
    intensity/=max(float(intensity.max()),1e-30)
    return {"method":"2D scalar TE finite-difference frequency domain; reciprocal waveguide-mode launch",
        "model":asdict(model),"grid":{"nx":x.size,"nz":z.size,"cells":x.size*z.size,"dx_um":dx,"dz_um":dz},
        "mode":{"n_eff":beta/k0,"beta_per_um":beta,"normalization":"integral |Ey|^2 dz = 1"},
        "efficiencies":{"target_free_space_mode":target_eff,"upward_radiation":upward,
            "substrate_radiation":substrate,"residual_forward_waveguide":residual_guided,
            "back_reflection":back_reflection,"accounted_power":accounted,
            "numerical_or_absorber_residual":1-accounted,
            "insertion_loss_db":float(-10*np.log10(max(target_eff,1e-30))),
            "directionality":float(upward/max(upward+substrate,1e-30))},
        "field":{"x_um":x[::stride].tolist(),"z_um":z[::stride].tolist(),
                 "normalized_Ey2":intensity.T.tolist()},
        "reciprocity_statement":"For reciprocal isotropic materials, power emitted by the launched waveguide mode into the normalized target beam equals coupling from the time-reversed target beam into that waveguide mode.",
        "scope":"Initial 2D TE finite-device solver. It includes finite length, partial etch, substrate leakage, back-reflection, absorbing boundaries, and output-mode normalization. It does not include finite lateral width or full-vector 3D polarization mixing.",
        "references":[{"title":"Optical Waveguide Theory","authors":"Snyder and Love","applies_to":"mode normalization and reciprocity"},
          {"title":"Grating couplers for coupling between optical fibers and nanophotonic waveguides","doi":"10.1143/JJAP.45.6071"},
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
    power_residuals={name:abs(run["efficiencies"]["numerical_or_absorber_residual"])
                     for name,run in (("coarse",coarse),("refined",refined),
                                      ("absorber",absorber),("padding",padding))}
    maximum=max([*changes.values(),*absorber_changes.values(),*padding_changes.values()])
    research=maximum <= 1e-2 and max(power_residuals.values()) <= 1e-2
    publication=maximum <= 2e-3 and max(power_residuals.values()) <= 2e-3
    return {"coarse":coarse,"refined":refined,"absorber_variant":absorber,
            "padding_variant":padding,"changes":changes,
            "absorber_changes":absorber_changes,"padding_changes":padding_changes,
            "power_budget_residuals":power_residuals,"maximum_change":maximum,
            "passed_research":research,"passed_publication_screen":publication,
            "interpretation":"The certificate varies mesh, absorber strength, and all domain paddings independently. Wavelength sampling and comparison with an external benchmark remain separate publication checks."}
