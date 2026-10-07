"""Two-branch Hermitian coupled-oscillator fit in photon-energy units."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

HC_EV_UM = 1.239841984


def fit_coupled_branches(parameter, branch_1_um, branch_2_um,
                         cavity_linewidth_mev: float | None = None,
                         matter_linewidth_mev: float | None = None) -> dict:
    x = np.asarray(parameter, dtype=float)
    l1 = np.asarray(branch_1_um, dtype=float)
    l2 = np.asarray(branch_2_um, dtype=float)
    if x.ndim != 1 or l1.shape != x.shape or l2.shape != x.shape or len(x) < 4:
        raise ValueError("Use at least four parameter values with two branch wavelengths each")
    if not np.isfinite([*x, *l1, *l2]).all() or np.any(l1 <= 0) or np.any(l2 <= 0):
        raise ValueError("Branch parameters and wavelengths must be positive/finite where applicable")
    if np.ptp(x) <= 0:
        raise ValueError("Branch parameter values must span a nonzero interval")
    energy_a, energy_b = HC_EV_UM/l1, HC_EV_UM/l2
    observed_upper = np.maximum(energy_a, energy_b)
    observed_lower = np.minimum(energy_a, energy_b)
    center = float(np.mean(x))
    scale = max(float(np.ptp(x)), 1e-12)
    u = (x-center)/scale
    all_energy = np.concatenate((observed_lower, observed_upper))
    span = max(float(np.ptp(all_energy)), 1e-4)
    exciton0 = float(np.mean((observed_upper+observed_lower)/2))
    cavity0 = exciton0
    slope0 = float(((observed_upper[-1]+observed_lower[-1])
                    -(observed_upper[0]+observed_lower[0]))/2)
    g0 = max(float(np.min(observed_upper-observed_lower))/2, span*.02)

    def predict(values):
        exciton, cavity_at_center, cavity_slope, coupling = values
        cavity = cavity_at_center+cavity_slope*u
        detuning = cavity-exciton
        root = np.sqrt((detuning/2)**2+coupling**2)
        return cavity, detuning, (cavity+exciton)/2+root, (cavity+exciton)/2-root

    def residual(values):
        _, _, upper, lower = predict(values)
        return np.concatenate((upper-observed_upper, lower-observed_lower))

    lo = float(np.min(all_energy)-2*span)
    hi = float(np.max(all_energy)+2*span)
    fit = least_squares(residual, [exciton0, cavity0, slope0, g0],
                        bounds=([lo, lo, -10*span, 0],
                                [hi, hi, 10*span, 5*span]),
                        max_nfev=20000, x_scale="jac")
    exciton, cavity_at_center, cavity_slope, coupling = fit.x
    cavity, detuning, predicted_upper, predicted_lower = predict(fit.x)
    denominator = np.sqrt(detuning**2+4*coupling**2)
    upper_photonic = .5*(1+detuning/np.maximum(denominator, 1e-300))
    lower_photonic = 1-upper_photonic
    resonance_u = (exciton-cavity_at_center)/cavity_slope if abs(cavity_slope) > 1e-14 else np.nan
    resonance_parameter = center+scale*resonance_u if np.isfinite(resonance_u) else None
    rms = float(np.sqrt(np.mean(residual(fit.x)**2)))
    result = {
        "parameter": x.tolist(), "parameter_center": center, "parameter_span": scale,
        "observed_upper_energy_ev": observed_upper.tolist(),
        "observed_lower_energy_ev": observed_lower.tolist(),
        "predicted_upper_energy_ev": predicted_upper.tolist(),
        "predicted_lower_energy_ev": predicted_lower.tolist(),
        "predicted_upper_wavelength_um": (HC_EV_UM/predicted_upper).tolist(),
        "predicted_lower_wavelength_um": (HC_EV_UM/predicted_lower).tolist(),
        "bare_cavity_energy_ev": cavity.tolist(),
        "bare_matter_energy_ev": float(exciton),
        "detuning_ev": detuning.tolist(),
        "coupling_mev": float(1000*coupling),
        "minimum_splitting_mev": float(2000*coupling),
        "resonance_parameter": None if resonance_parameter is None else float(resonance_parameter),
        "upper_photonic_fraction": upper_photonic.tolist(),
        "upper_matter_fraction": (1-upper_photonic).tolist(),
        "lower_photonic_fraction": lower_photonic.tolist(),
        "lower_matter_fraction": (1-lower_photonic).tolist(),
        "rmse_mev": 1000*rms, "success": bool(fit.success),
        "model": "Hermitian 2×2 coupled oscillator with one linear bare photonic branch and one constant bare matter energy",
        "interpretation": "Hopfield fractions describe fitted eigenvector composition. They are not illumination, collection, or power-transfer efficiencies."
    }
    if cavity_linewidth_mev is not None or matter_linewidth_mev is not None:
        if cavity_linewidth_mev is None or matter_linewidth_mev is None:
            raise ValueError("Provide both bare cavity and matter linewidths, or leave both blank")
        gc, gx = float(cavity_linewidth_mev), float(matter_linewidth_mev)
        if not np.isfinite([gc, gx]).all() or min(gc, gx) <= 0:
            raise ValueError("Bare linewidths must be positive finite values in meV")
        result["linewidths_mev"] = {"cavity": gc, "matter": gx}
        result["cooperativity_4g2_over_product"] = float(4*(1000*coupling)**2/(gc*gx))
        result["splitting_to_mean_linewidth"] = float((2000*coupling)/((gc+gx)/2))
    return result
