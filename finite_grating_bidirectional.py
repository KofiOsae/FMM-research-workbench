"""Independent Gaussian-incident solve and bidirectional reciprocity checks.

This module shares the physical grid, scalar TE operator, material raster and
guided-mode projection used by :mod:`finite_grating`.  It deliberately uses a
different source construction from the guided-port launch: a downward
Gaussian angular spectrum is converted to an exact discrete equivalent
current through a homogeneous-cladding reference operator.
"""
from __future__ import annotations

from dataclasses import asdict, replace
import numpy as np
from scipy.sparse.linalg import spsolve


def _complex_value(value: complex) -> dict:
    return {"real": float(np.real(value)), "imag": float(np.imag(value)),
            "magnitude": float(abs(value)),
            "phase_deg": float(np.rad2deg(np.angle(value)))}


def _target_line(model, x, k0):
    from finite_grating import _grating_length
    center = (model.target_center_um if model.target_center_um is not None
              else _grating_length(model)/2)
    kx0 = model.cladding_n*k0*np.sin(np.deg2rad(model.target_angle_deg))
    outgoing = np.exp(-((x-center)/model.target_waist_um)**2)*np.exp(
        1j*(kx0*(x-center)+np.deg2rad(model.target_phase_deg)))
    return center, kx0, outgoing


def _field_payload(total, incident, scattered, eps, x, z, k0):
    stride=max(1,int(max(x.size,z.size)/240)); sampled=total[::stride,::stride]
    scale=max(float(np.max(abs(sampled))),1e-30)
    dx=float(x[1]-x[0]); dz=float(z[1]-z[0])
    d_ey_dx=np.gradient(total,dx,axis=0); d_ey_dz=np.gradient(total,dz,axis=1)
    hx=1j*d_ey_dz/k0; hz=-1j*d_ey_dx/k0
    sx=.5*np.real(total*np.conj(hz)); sz=-.5*np.real(total*np.conj(hx))
    pscale=max(float(np.max(np.sqrt(sx*sx+sz*sz))),1e-30)
    return {"x_um":x[::stride].tolist(),"z_um":z[::stride].tolist(),
        "normalization":"Every electric-field view uses the total-field max |Ey|; relative Poynting components use max sqrt(Sx^2+Sz^2).",
        "normalized_Ey2":(abs(sampled/scale)**2).T.tolist(),
        "normalized_Ey_real":(sampled.real/scale).T.tolist(),
        "normalized_Ey_imag":(sampled.imag/scale).T.tolist(),
        "Ey_phase_rad":np.angle(sampled).T.tolist(),
        "incident_Ey2":(abs(incident[::stride,::stride]/scale)**2).T.tolist(),
        "incident_Ey_real":(incident[::stride,::stride].real/scale).T.tolist(),
        "scattered_Ey2":(abs(scattered[::stride,::stride]/scale)**2).T.tolist(),
        "scattered_Ey_real":(scattered[::stride,::stride].real/scale).T.tolist(),
        "normalized_Sx":(sx[::stride,::stride]/pscale).T.tolist(),
        "normalized_Sz":(sz[::stride,::stride]/pscale).T.tolist(),
        "epsilon_r":eps[::stride,::stride].real.T.tolist()}


def solve_fiber_incident(model):
    """Solve Gaussian port -> finite grating -> guided waveguide mode.

    Angle is expressed in global coordinates.  The reciprocal incident port is
    the time reverse of the requested outgoing Gaussian, so its transverse
    wavevector and signed global angle have the opposite sign.
    """
    from finite_grating import (_axes, _permittivity, _absorber_sigma,
        _operator, _te_mode, _mode_amplitudes, _angular_power,
        _angular_spectrum_rows, _raster_geometry_metrics)

    model.validate(); x,z=_axes(model); dx=float(x[1]-x[0]); dz=float(z[1]-z[0])
    eps,reference=_permittivity(model,x,z); k0=2*np.pi/model.wavelength_um
    beta,phi=_te_mode(reference,z,k0)
    sigma_x=_absorber_sigma(x,x[0],x[-1],model.absorber_um,model.absorber_strength)
    sigma_z=_absorber_sigma(z,z[0],z[-1],model.absorber_um,model.absorber_strength)
    sx=1+1j*sigma_x; sz=1+1j*sigma_z
    actual_op=_operator(eps,dx,dz,sx,sz,k0)
    background_eps=np.full(eps.shape,model.cladding_n**2,dtype=complex)
    background_op=_operator(background_eps,dx,dz,sx,sz,k0)

    left_i=max(2,int(round(model.absorber_um/dx))+3)
    source_i=min(x.size-4,int(round((x[-1]-model.absorber_um-2*dx-x[0])/dx)))
    right_i=max(left_i+3,source_i-4)
    top_j=min(z.size-5,int(round((z[-1]-model.absorber_um-5*dz-z[0])/dz)))
    bottom_j=max(2,int(round(model.absorber_um/dz))+3)

    _,_,outgoing=_target_line(model,x,k0)
    # Lorentz-reciprocal incoming port: complex conjugate of the outgoing
    # transverse field and negative longitudinal wavevector.
    line=np.conj(outgoing)
    spectrum=dx*np.fft.fftshift(np.fft.fft(np.fft.ifftshift(line)))
    kx=2*np.pi*np.fft.fftshift(np.fft.fftfreq(x.size,d=dx))
    kz2=(model.cladding_n*k0)**2-kx**2; prop=kz2>0
    kz=np.sqrt(np.maximum(kz2,0)); spectrum[~prop]=0
    injection_j=min(z.size-3,top_j+3); z_ref=z[top_j]
    desired=np.zeros(eps.shape,complex)
    for j in range(1,injection_j+1):
        propagated=spectrum*np.exp(-1j*kz*(z[j]-z_ref))
        desired[:,j]=np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(propagated/dx)))

    # Truncate only inside the absorbing rim.  Applying the discrete
    # background operator produces the equivalent current, including the
    # electric/magnetic Huygens-sheet pair implicit in the stencil.
    wx=np.ones(x.size); wz=np.ones(z.size)
    ex=max(2,int(np.ceil(model.absorber_um/dx)))
    ez=max(2,int(np.ceil(model.absorber_um/dz)))
    rx=np.sin(np.linspace(0,np.pi/2,ex))**2
    rz=np.sin(np.linspace(0,np.pi/2,ez))**2
    wx[:ex]=rx; wx[-ex:]=rx[::-1]; wz[:ez]=rz; wz[-ez:]=rz[::-1]
    desired*=wx[:,None]*wz[None,:]
    desired[[0,-1],:]=0; desired[:,[0,-1]]=0
    rhs=background_op@desired.ravel()
    incident=spsolve(background_op,rhs).reshape(eps.shape)
    total=spsolve(actual_op,rhs).reshape(eps.shape)
    scattered=total-incident

    *_,incident_power,_=_angular_power(incident,top_j,x,z,model.cladding_n,k0,False)
    if incident_power <= 1e-20:
        raise ValueError("The incident Gaussian has negligible propagating power; increase the waist or domain width")
    leftward,_=_mode_amplitudes(total,left_i,phi,z,dx,beta)
    _,rightward=_mode_amplitudes(total,right_i,phi,z,dx,beta)
    guided_left=beta*abs(leftward)**2/incident_power
    guided_right=beta*abs(rightward)**2/incident_power
    coefficient=np.sqrt(beta)*rightward/np.sqrt(incident_power)
    rkx,rkz,rup,rprop,reflected_power,rdk=_angular_power(
        scattered,top_j,x,z,model.cladding_n,k0,True)
    lower_n=(model.handle_n if model.box_thickness_um is not None and model.handle_n is not None
             else model.substrate_n)
    *_,substrate_power,_=_angular_power(total,bottom_j,x,z,lower_n,k0,False)
    reflected=reflected_power/incident_power; substrate=substrate_power/incident_power
    accounted=guided_left+guided_right+reflected+substrate

    akx,akz,down,aprop,_,_= _angular_power(incident,top_j,x,z,model.cladding_n,k0,False)
    weights=akz[aprop]*abs(down[aprop])**2
    mean_kx=float(np.sum(akx[aprop]*weights)/max(np.sum(weights),1e-30))
    actual_angle=float(np.rad2deg(np.arcsin(np.clip(mean_kx/(model.cladding_n*k0),-1,1))))
    e={"selected_mode_coupling":float(guided_right),
       "target_free_space_mode":float(guided_right), # compatibility alias
       "guided_output_mode":float(guided_right),
       "upward_radiation":float(reflected),"substrate_radiation":float(substrate),
       "residual_forward_waveguide":float(guided_left),"back_reflection":float(reflected),
       "accounted_power":float(accounted),"numerical_or_absorber_residual":float(1-accounted),
       "insertion_loss_db":float(-10*np.log10(max(guided_right,1e-30))),
       "directionality":float(guided_right/max(guided_left+guided_right,1e-30))}
    core=(z>=0)&(z<=model.waveguide_height_um)
    field=_field_payload(total,incident,scattered,eps,x,z,k0)
    field["annotations"]={"source":{"axis":"z","position_um":float(z[injection_j]),
        "label":"Gaussian equivalent-current plane"},
        "monitors":[{"axis":"x","position_um":float(x[left_i]),"label":"left guided port"},
                    {"axis":"x","position_um":float(x[right_i]),"label":"right guided port"},
                    {"axis":"z","position_um":float(z[top_j]),"label":"incident/reflection monitor"},
                    {"axis":"z","position_um":float(z[bottom_j]),"label":"lower radiation monitor"}],
        "absorber_um":float(model.absorber_um)}
    stack=({"type":"semi_infinite_substrate","substrate_n":model.substrate_n}
           if model.box_thickness_um is None else
           {"type":"finite_box_and_handle","box_n":model.substrate_n,
            "box_thickness_um":model.box_thickness_um,"handle_n":model.handle_n})
    return {"method":"2D scalar TE FDFD; independent Gaussian-incident in-coupling solve",
      "model":asdict(model),"grid":{"nx":x.size,"nz":z.size,"cells":x.size*z.size,"dx_um":dx,"dz_um":dz},
      "geometry":{"upper_cladding_n":model.cladding_n,"device_n":model.core_n,
          "lower_stack":stack,"rasterization":model.discretization,
          "subpixel_samples_per_axis":model.subpixel_samples if model.discretization=="cell_average" else 1,
          "nominal_and_rasterized":_raster_geometry_metrics(model,x,z,eps,reference)},
      "formulation":{"time_convention":"exp(-i omega t)",
          "equation":"(d_x^2 + d_z^2 + k0^2 epsilon_r) Ey = b",
          "source":"discrete equivalent current b=A_cladding E_inc for a downward Gaussian angular spectrum",
          "excitation":"gaussian_port_to_grating_to_waveguide",
          "incident_power_relative":float(incident_power),
          "requested_outgoing_angle_deg":float(model.target_angle_deg),
          "reciprocal_incoming_global_angle_deg":float(-model.target_angle_deg),
          "sampled_incoming_global_angle_deg":actual_angle,
          "injection_z_um":float(z[injection_j]),"monitor_z_um":float(z[top_j]),
          "reference_reproduction_relative_error":float(np.linalg.norm(incident-desired)/max(np.linalg.norm(desired),1e-30))},
      "mode":{"n_eff":beta/k0,"beta_per_um":beta,
          "device_core_fraction":float(np.trapezoid(abs(phi[core])**2,z[core])),
          "selection":"maximum electric-profile localization in the device layer above z=0",
          "normalization":"integral |Ey|^2 dz = 1"},
      "efficiencies":e,"selected_port_coefficient":_complex_value(coefficient),
      "angular_spectrum":{"global_angle_convention":"Outgoing angle is measured from global +z toward +x; the reciprocal incoming wave reverses the full wavevector.",
          "incident":_angular_spectrum_rows(akx,akz,down,aprop,model.cladding_n,k0,incident_power),
          "reflected":_angular_spectrum_rows(rkx,rkz,rup,rprop,model.cladding_n,k0,incident_power),
          "density_integral_note":"Integrate power_fraction_density_per_rad_per_um over kx (rad/µm) to recover the channel fraction."},
      "field":field,
      "reciprocity_statement":"This is a genuine, independent Gaussian-source field solve. Compare it with the separate guided-source solve using the bidirectional reciprocity certificate.",
      "scope":"2D scalar TE finite-device in-coupling with a Gaussian angular-spectrum port. The beam is a width-invariant Gaussian sheet, not a full 3D circular fiber mode. Full-vector polarization mixing, finite lateral width and multiple guided modes are outside this model.",
      "references":[{"title":"Optical Waveguide Theory","authors":"Snyder and Love","applies_to":"mode normalization and reciprocity"},
        {"title":"Grating couplers for coupling between optical fibers and nanophotonic waveguides","doi":"10.1143/JJAP.45.6071"},
        {"title":"Computational Electrodynamics: The Finite-Difference Time-Domain Method","authors":"Taflove and Hagness","applies_to":"equivalent-current and total/scattered-field source construction"}]}


def reciprocity_certificate(model, refine: bool = False):
    """Run two independent sources and compare the paired port coefficients.

    When ``refine`` is true, repeat both directions at a 1.25-times finer mesh
    and require the two directional efficiencies themselves to stabilize.
    """
    from finite_grating import solve_finite_grating
    forward=solve_finite_grating(replace(model,excitation="waveguide"))
    reverse=solve_fiber_incident(replace(model,excitation="fiber"))
    f=forward["selected_port_coefficient"]; r=reverse["selected_port_coefficient"]
    fc=complex(f["real"],f["imag"]); rc=complex(r["real"],r["imag"])
    # Reference planes differ by a known but currently uncalibrated propagation
    # phase, so magnitude/efficiency agreement is the primary invariant.  The
    # raw complex coefficients remain exported for phase-reference auditing.
    magnitude_error=abs(abs(fc)-abs(rc))/max(abs(fc),abs(rc),1e-30)
    efficiency_error=abs(abs(fc)**2-abs(rc)**2)
    report={"forward":forward,"reverse":reverse,
        "forward_complex_coefficient":f,"reverse_complex_coefficient":r,
        "relative_magnitude_error":float(magnitude_error),
        "absolute_efficiency_error":float(efficiency_error),
        "passed_research_screen":bool(magnitude_error<=.05 and efficiency_error<=.02),
        "criterion":{"maximum_relative_magnitude_error":.05,
                     "maximum_absolute_efficiency_error":.02},
        "phase_status":"Raw phases are reported, but phase reciprocity is not certified until both port reference planes are de-embedded to the same origin.",
        "interpretation":"Two independent linear solves are compared: guided eigenmode launch and a time-reversed Gaussian equivalent-current launch. A failed screen can reflect discretization, absorber, finite-window port truncation, or inconsistent reference planes; it must not be hidden by averaging the two directions."}
    if refine:
        refined_model=replace(model,mesh_um=model.mesh_um/1.25)
        refined=reciprocity_certificate(refined_model,False)
        forward_change=abs(refined["forward"]["efficiencies"]["selected_mode_coupling"]-
                           forward["efficiencies"]["selected_mode_coupling"])
        reverse_change=abs(refined["reverse"]["efficiencies"]["selected_mode_coupling"]-
                           reverse["efficiencies"]["selected_mode_coupling"])
        refinement_pass=bool(max(forward_change,reverse_change)<=.02 and
                             refined["relative_magnitude_error"]<=.05 and
                             refined["absolute_efficiency_error"]<=.02)
        report["mesh_refinement"]={"factor":1.25,"mesh_um":refined_model.mesh_um,
            "forward_efficiency_change":float(forward_change),
            "reverse_efficiency_change":float(reverse_change),
            "refined_relative_magnitude_error":refined["relative_magnitude_error"],
            "refined_absolute_efficiency_error":refined["absolute_efficiency_error"],
            "passed":refinement_pass,
            "criterion":"Both directional efficiencies change by at most 0.02 absolute and the refined reciprocity errors pass the base screen."}
        report["refined_coefficients"]={"forward":refined["forward_complex_coefficient"],
                                        "reverse":refined["reverse_complex_coefficient"]}
        report["passed_research_screen"]=bool(report["passed_research_screen"] and refinement_pass)
    return report
