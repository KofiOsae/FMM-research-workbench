"""Bound modes of an isotropic, planar three-medium dielectric waveguide.

The guided wave propagates along +x. The core occupies 0 <= z <= thickness.
This is a separate eigenmode calculation from normally incident RCWA.
"""

from dataclasses import asdict, dataclass, replace

import numpy as np
from scipy.optimize import brentq

from materials import material_n


@dataclass(frozen=True)
class WaveguideModel:
    wavelength_um: float = 1.55
    thickness_um: float = .35
    core_material: str = "dielectric"
    core_n: float = 2.0
    top_material: str = "air"
    top_n: float = 1.0
    bottom_material: str = "silica"
    bottom_n: float = 1.45
    polarization: str = "TE"

    def indices(self) -> tuple[float, float, float]:
        if not np.isfinite([self.wavelength_um, self.thickness_um]).all() or self.wavelength_um <= 0 or self.thickness_um <= 0:
            raise ValueError("Waveguide wavelength and core thickness must be finite and positive")
        result = []
        for material, constant in ((self.core_material, self.core_n),
                                   (self.top_material, self.top_n),
                                   (self.bottom_material, self.bottom_n)):
            n = material_n(material, self.wavelength_um, constant)
            if abs(n.imag) > 1e-9:
                raise ValueError("This guided-mode solver accepts lossless isotropic dielectrics only")
            result.append(float(n.real))
        core, top, bottom = result
        if core <= max(top, bottom):
            raise ValueError("A bound mode needs a core index above both cladding indices")
        if self.polarization not in ("TE", "TM"):
            raise ValueError("Guided-mode polarization must be TE or TM")
        return core, top, bottom


def _roots(model: WaveguideModel):
    core, top, bottom = model.indices()
    k0 = 2*np.pi/model.wavelength_um
    lower, upper = max(top, bottom), core

    def terms(neff):
        h = k0*np.sqrt(max(0, core**2-neff**2))
        qt = k0*np.sqrt(max(0, neff**2-top**2))
        qb = k0*np.sqrt(max(0, neff**2-bottom**2))
        wt = core**2/top**2 if model.polarization == "TM" else 1
        wb = core**2/bottom**2 if model.polarization == "TM" else 1
        return h, qt, qb, np.arctan2(wt*qt, h), np.arctan2(wb*qb, h)

    def equation(neff, mode):
        h, _, _, pt, pb = terms(neff)
        return h*model.thickness_um-pt-pb-mode*np.pi

    eps = min((upper-lower)*1e-9, 1e-10)
    roots = []
    for mode in range(50):
        if equation(lower+eps, mode) <= 0:
            break
        if equation(upper-eps, mode) >= 0:
            raise RuntimeError("Guided-mode root was not bracketed")
        neff = brentq(equation, lower+eps, upper-eps, args=(mode,), xtol=1e-13)
        roots.append((mode, neff, terms(neff)))
    if equation(lower+eps, 50) > 0:
        raise ValueError("More than 50 guided modes are present; reduce the core thickness for this interactive solver")
    return roots, (core, top, bottom)


def solve_waveguide(model: WaveguideModel) -> dict:
    roots, (core, top, bottom) = _roots(model)
    modes = []
    for order, neff, (h, qt, qb, pt, pb) in roots:
        extent_b = min(12/max(qb,1e-12), 100*model.thickness_um)
        extent_t = min(12/max(qt,1e-12), 100*model.thickness_um)
        z = np.linspace(-extent_b, model.thickness_um+extent_t, 1201)
        profile = np.where(z < 0, np.cos(pb)*np.exp(qb*np.minimum(z,0)),
                           np.where(z > model.thickness_um,
                                    np.cos(h*model.thickness_um-pb)*np.exp(-qt*np.maximum(z-model.thickness_um,0)),
                                    np.cos(h*z-pb)))
        intensity = profile**2
        if model.polarization == "TM":
            # TM plotted quantity is magnetic H_y intensity, not electric intensity.
            component = "|H_y|²"
        else:
            component = "|E_y|²"
        # Exact integral of the squared scalar mode over the infinite claddings.
        core_integral = (model.thickness_um/2
                         + (np.sin(2*(h*model.thickness_um-pb))+np.sin(2*pb))/(4*h))
        bottom_integral = np.cos(pb)**2/(2*qb)
        top_integral = np.cos(h*model.thickness_um-pb)**2/(2*qt)
        profile_fraction = core_integral/(core_integral+bottom_integral+top_integral)
        if model.polarization == "TM":
            # Along-waveguide Poynting flux is proportional to |H_y|²/epsilon.
            guided_fraction = (core_integral/core**2
                               / (core_integral/core**2 + bottom_integral/bottom**2
                                  + top_integral/top**2))
        else:
            # For TE, along-waveguide flux is proportional to |E_y|².
            guided_fraction = profile_fraction
        group_index = None
        delta = model.wavelength_um*1e-4
        group_index_change = None
        try:
            minus = _roots(replace(model,wavelength_um=model.wavelength_um-delta))[0][order][1]
            plus = _roots(replace(model,wavelength_um=model.wavelength_um+delta))[0][order][1]
            group_index = float(neff-model.wavelength_um*(plus-minus)/(2*delta))
            half = delta/2
            minus_half = _roots(replace(model,wavelength_um=model.wavelength_um-half))[0][order][1]
            plus_half = _roots(replace(model,wavelength_um=model.wavelength_um+half))[0][order][1]
            group_index_half = neff-model.wavelength_um*(plus_half-minus_half)/(2*half)
            group_index_change = float(abs(group_index-group_index_half))
        except (ValueError, IndexError):
            pass
        modes.append({"order":order,"name":f"{model.polarization}{order}",
                      "n_eff":float(neff),"beta_rad_um":float(2*np.pi*neff/model.wavelength_um),
                      "group_index":group_index,"group_index_step_change":group_index_change,
                      "confinement_fraction":float(guided_fraction),
                      "profile_integral_fraction":float(profile_fraction),
                      "top_decay_um":float(1/qt),"bottom_decay_um":float(1/qb),
                      "boundary_residual":float(abs(h*model.thickness_um-pt-pb-order*np.pi)),
                      "z_um":z.tolist(),"profile":(intensity/intensity.max()).tolist(),
                      "profile_quantity":component})
    return {"model":asdict(model),"indices":{"core":core,"top":top,"bottom":bottom},
            "modes":modes,"note":"Guided modes propagate in-plane; this solver uses lossless planar media and does not predict coupling efficiency from the incident plane wave."}
