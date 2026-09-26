"""Electric-dipole decay rate above a laterally uniform multilayer.

Implements the Chance--Prock--Silbey / dyadic Green-function Sommerfeld
integrals.  The emitter is in the top semi-infinite, real-index medium and the
finite layers are ordered from top to bottom as in :mod:`stack`.
"""

from __future__ import annotations

from dataclasses import asdict, replace

import numpy as np
from scipy.integrate import quad_vec

from materials import material_n
from stack import StackModel


def _forward_sqrt(value: complex | np.ndarray) -> complex | np.ndarray:
    root = np.lib.scimath.sqrt(value)
    return np.where((np.imag(root) < 0) |
                    ((np.abs(np.imag(root)) < 1e-14) & (np.real(root) < 0)),
                    -root, root)


def multilayer_reflection(model: StackModel, u: np.ndarray,
                          polarization: str) -> np.ndarray:
    """Reflection seen by a source in the top medium versus k_parallel/k_top.

    The p coefficient follows the Green-function convention (r_p -> +1 for a
    perfect electric conductor); its normal-incidence sign is opposite to the
    tangential-E Fresnel coefficient used by some transfer-matrix texts.
    """
    if polarization not in ("s", "p"):
        raise ValueError("Polarization must be s or p")
    wavelength = model.wavelength_um
    indices = [complex(model.incident_n)]
    indices += [material_n(layer.background_material, wavelength,
                           layer.background_n) for layer in model.layers]
    indices += [complex(model.exit_n)]
    eps = np.asarray(indices, dtype=complex) ** 2
    k0 = 2*np.pi/wavelength
    k_parallel = k0*indices[0]*np.asarray(u, dtype=complex)
    kz = [_forward_sqrt((k0*n)**2-k_parallel**2) for n in indices]
    reflection = np.zeros_like(k_parallel, dtype=complex)
    for j in range(len(indices)-2, -1, -1):
        if polarization == "s":
            interface = (kz[j]-kz[j+1])/(kz[j]+kz[j+1])
        else:
            interface = (eps[j+1]*kz[j]-eps[j]*kz[j+1]) / \
                        (eps[j+1]*kz[j]+eps[j]*kz[j+1])
        phase = (np.exp(2j*kz[j+1]*model.layers[j].thickness_um)
                 if j < len(model.layers) else 1+0j)
        reflection = (interface+reflection*phase)/(1+interface*reflection*phase)
    return reflection


def _compute(model: StackModel, distance_um: float, points: int,
             evanescent_limit: float) -> dict:
    # Variable changes remove the light-line square-root singularity. Adaptive
    # vector quadrature then resolves narrow surface-mode poles.
    theta = np.linspace(1e-7, np.pi/2-1e-7, points)
    u_prop = np.sin(theta)
    rp = multilayer_reflection(model, u_prop, "p")
    rs = multilayer_reflection(model, u_prop, "s")
    k = 2*np.pi*model.incident_n/model.wavelength_um
    phase = np.exp(2j*k*distance_um*np.cos(theta))
    perp_prop_density = np.real(np.sin(theta)**3*rp*phase)
    para_prop_density = np.real(np.sin(theta)*(rs-np.cos(theta)**2*rp)*phase)
    def prop_integrand(angle):
        uu = np.asarray([np.sin(angle)])
        rpv = multilayer_reflection(model, uu, "p")[0]
        rsv = multilayer_reflection(model, uu, "s")[0]
        ph = np.exp(2j*k*distance_um*np.cos(angle))
        return np.asarray([1.5*np.real(np.sin(angle)**3*rpv*ph),
                           .75*np.real(np.sin(angle)*(rsv-np.cos(angle)**2*rpv)*ph)])
    prop, _ = quad_vec(prop_integrand, 0, np.pi/2, epsabs=2e-7,
                       epsrel=2e-6, limit=max(100, points))

    t = np.linspace(1e-7, evanescent_limit, points)
    u_evan = np.sqrt(1+t*t)
    rp_e = multilayer_reflection(model, u_evan, "p")
    rs_e = multilayer_reflection(model, u_evan, "s")
    decay = np.exp(-2*k*distance_um*t)
    perp_evan_density = np.real(-1j*u_evan**2*rp_e*decay)
    para_evan_density = np.real(-1j*(rs_e+t*t*rp_e)*decay)
    def evan_integrand(tt):
        uu = np.asarray([np.sqrt(1+tt*tt)])
        rpv = multilayer_reflection(model, uu, "p")[0]
        rsv = multilayer_reflection(model, uu, "s")[0]
        dec = np.exp(-2*k*distance_um*tt)
        return np.asarray([1.5*np.real(-1j*uu[0]**2*rpv*dec),
                           .75*np.real(-1j*(rsv+tt*tt*rpv)*dec)])
    evan, _ = quad_vec(evan_integrand, 0, evanescent_limit, epsabs=2e-7,
                       epsrel=2e-6, limit=max(100, points))
    perp_prop, para_prop = prop
    perp_evan, para_evan = evan
    return {
        "perpendicular": 1+float(perp_prop+perp_evan),
        "parallel": 1+float(para_prop+para_evan),
        "perpendicular_propagating_correction": float(perp_prop),
        "perpendicular_evanescent_correction": float(perp_evan),
        "parallel_propagating_correction": float(para_prop),
        "parallel_evanescent_correction": float(para_evan),
        "u_propagating": u_prop.tolist(),
        "perpendicular_propagating_density": perp_prop_density.tolist(),
        "parallel_propagating_density": para_prop_density.tolist(),
        "u_evanescent": u_evan.tolist(),
        "perpendicular_evanescent_density": perp_evan_density.tolist(),
        "parallel_evanescent_density": para_evan_density.tolist(),
    }


def dipole_ldos(model: StackModel, distance_um: float, orientation: str = "isotropic",
                points: int = 240, evanescent_limit: float = 30.0,
                collection_na: float | None = None) -> dict:
    """Return normalized electric LDOS / decay rate for a dipole above a stack."""
    model.validate()
    if any(layer.kind != "uniform" for layer in model.layers):
        raise ValueError("Dipole LDOS currently requires uniform layers; patterned Green functions need a periodic-source solver")
    if abs(complex(model.incident_n).imag) > 1e-14 or model.incident_n <= 0:
        raise ValueError("The emitter medium must have a positive real refractive index")
    if not np.isfinite([distance_um, evanescent_limit]).all() or distance_um <= 0:
        raise ValueError("Emitter distance must be positive and finite")
    if orientation not in ("perpendicular", "parallel", "isotropic"):
        raise ValueError("Orientation must be perpendicular, parallel, or isotropic")
    if not 40 <= points <= 800 or not 2 <= evanescent_limit <= 200:
        raise ValueError("Use 40–800 quadrature points and an evanescent limit from 2 to 200")
    if collection_na is None:
        collection_na = float(model.incident_n)
    if not np.isfinite(collection_na) or not 0 < collection_na <= model.incident_n:
        raise ValueError("Collection NA must be positive and no larger than the emitter-medium index")
    result = _compute(model, distance_um, points, evanescent_limit)
    refined = _compute(model, distance_um, min(800, 2*points),
                       min(200, 1.25*evanescent_limit))
    result["isotropic"] = (result["perpendicular"]+2*result["parallel"])/3
    refined["isotropic"] = (refined["perpendicular"]+2*refined["parallel"])/3
    selected = result[orientation]
    change = abs(refined[orientation]-selected)
    tolerance = max(1e-3, 2e-3*abs(refined[orientation]))
    theta_max = float(np.arcsin(collection_na/model.incident_n))
    def collected(angle):
        u = np.asarray([np.sin(angle)])
        rp = multilayer_reflection(model, u, "p")[0]
        rs = multilayer_reflection(model, u, "s")[0]
        w = np.cos(angle); phase = np.exp(2j*(2*np.pi*model.incident_n/model.wavelength_um)*distance_um*w)
        perpendicular = .75*np.sin(angle)**3*abs(1+rp*phase)**2
        parallel = .375*np.sin(angle)*(abs(1+rs*phase)**2+w*w*abs(1-rp*phase)**2)
        return np.asarray([perpendicular, parallel], dtype=float)
    collected_values, _ = quad_vec(collected, 0, theta_max, epsabs=2e-7, epsrel=2e-6)
    c_perp, c_para = map(float, collected_values)
    c_iso = (c_perp+2*c_para)/3
    c_selected = {"perpendicular": c_perp, "parallel": c_para, "isotropic": c_iso}[orientation]
    return {
        **result,
        "orientation": orientation,
        "normalized_decay_rate": selected,
        "convergence": {"refined_value": refined[orientation],
                        "absolute_change": change, "tolerance": tolerance,
                        "converged": bool(change <= tolerance),
                        "coarse_points": points,
                        "refined_points": min(800, 2*points),
                        "coarse_evanescent_limit": evanescent_limit,
                        "refined_evanescent_limit": min(200, 1.25*evanescent_limit)},
        "distance_um": distance_um,
        "collection": {"NA": collection_na, "half_angle_deg": float(np.degrees(theta_max)),
            "perpendicular_power_over_homogeneous_total": c_perp,
            "parallel_power_over_homogeneous_total": c_para,
            "isotropic_power_over_homogeneous_total": c_iso,
            "selected_collection_efficiency": c_selected/max(selected, 1e-30),
            "definition": "power emitted into the upper-medium cone divided by total dipole power"},
        "wavelength_um": model.wavelength_um,
        "model": asdict(model),
        "method": "planar electric-dipole dyadic Green-function Sommerfeld integral",
        "normalization": "decay rate in the layered environment divided by the rate in the homogeneous emitter medium",
        "note": "The reported collection efficiency is rigorous for the selected upper-medium numerical aperture. Propagating and evanescent LDOS corrections are not by themselves a radiative/nonradiative split. A beta factor for a named guided or cavity mode requires projection onto that normalized mode."
    }


def dipole_ldos_spectrum(model: StackModel, start_um: float, stop_um: float,
                         wavelength_points: int, distance_um: float,
                         orientation: str = "isotropic", points: int = 120,
                         evanescent_limit: float = 30.0,
                         collection_na: float | None = None) -> dict:
    """Sweep planar LDOS using an independent, refined integral per wavelength."""
    if not np.isfinite([start_um, stop_um]).all() or start_um <= 0 or stop_um <= start_um:
        raise ValueError("LDOS spectrum requires 0 < start wavelength < stop wavelength")
    if not 3 <= wavelength_points <= 61:
        raise ValueError("Use 3–61 wavelengths for an LDOS spectrum")
    wavelengths = np.linspace(start_um, stop_um, wavelength_points)
    rows = []
    for wavelength in wavelengths:
        result = dipole_ldos(replace(model, wavelength_um=float(wavelength)),
                             distance_um, orientation, points,
                             evanescent_limit, collection_na)
        collection_key = f"{orientation}_power_over_homogeneous_total"
        rows.append({
            "wavelength_um": float(wavelength),
            "normalized_decay_rate": result["normalized_decay_rate"],
            "perpendicular": result["perpendicular"],
            "parallel": result["parallel"],
            "isotropic": result["isotropic"],
            "collection_efficiency": result["collection"]["selected_collection_efficiency"],
            "collection_power_over_homogeneous_total": result["collection"][collection_key],
            "converged": result["convergence"]["converged"],
            "refinement_change": result["convergence"]["absolute_change"],
        })
    return {
        "wavelength_um": wavelengths.tolist(), "rows": rows,
        "orientation": orientation, "distance_um": distance_um,
        "collection_na": collection_na if collection_na is not None else float(model.incident_n),
        "all_converged": all(row["converged"] for row in rows),
        "method": "independent planar dyadic Green-function Sommerfeld integration at every wavelength",
        "normalization": "decay rate in the layered environment divided by the homogeneous-emitter rate",
        "warning": "A total-LDOS peak can be lossy or evanescent. Collection is integrated separately over the upper-medium objective cone.",
    }
