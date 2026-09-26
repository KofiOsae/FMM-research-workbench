"""Local constitutive response and homogeneous anisotropic plane-wave modes.

This module deliberately does not pretend that a scalar RCWA stack accepts
tensor media.  It prepares and checks the tensor that a future vector solver
will consume and solves the exact homogeneous Maxwell eigenproblem.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eig


def _rotation_from_axis(theta_deg: float, phi_deg: float) -> np.ndarray:
    theta, phi = np.deg2rad([theta_deg, phi_deg])
    axis = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    trial = np.array([0., 0., 1.]) if abs(axis[2]) < .9 else np.array([1., 0., 0.])
    first = np.cross(trial, axis); first /= np.linalg.norm(first)
    second = np.cross(axis, first)
    return np.column_stack((first, second, axis))


def constitutive_response(wavelength_um: float, n_ordinary: float, k_ordinary: float,
                          n_extraordinary: float, k_extraordinary: float,
                          optic_axis_theta_deg: float = 0., optic_axis_phi_deg: float = 0.,
                          temperature_K: float = 293.15, reference_temperature_K: float = 293.15,
                          dn_dT_per_K: float = 0., intensity_W_m2: float = 0.,
                          n2_m2_W: float = 0., gain_per_m: float = 0.,
                          mu_relative: float = 1., propagation_theta_deg: float = 0.,
                          propagation_phi_deg: float = 0.) -> dict:
    values = np.asarray([wavelength_um, n_ordinary, k_ordinary, n_extraordinary,
        k_extraordinary, optic_axis_theta_deg, optic_axis_phi_deg, temperature_K,
        reference_temperature_K, dn_dT_per_K, intensity_W_m2, n2_m2_W,
        gain_per_m, mu_relative, propagation_theta_deg, propagation_phi_deg], float)
    if not np.isfinite(values).all() or wavelength_um <= 0 or min(n_ordinary, n_extraordinary, mu_relative) <= 0:
        raise ValueError("Use finite values and positive wavelength, real indices, and relative permeability")
    if temperature_K <= 0 or reference_temperature_K <= 0 or intensity_W_m2 < 0:
        raise ValueError("Temperatures must be positive and intensity cannot be negative")
    delta_n = dn_dT_per_K*(temperature_K-reference_temperature_K)+n2_m2_W*intensity_W_m2
    # Intensity gain g obeys I(z)=I(0)exp(gz), hence k_gain=-g*lambda/(4pi).
    gain_k = -gain_per_m*(wavelength_um*1e-6)/(4*np.pi)
    no = complex(n_ordinary+delta_n, k_ordinary+gain_k)
    ne = complex(n_extraordinary+delta_n, k_extraordinary+gain_k)
    principal = np.diag([no*no, no*no, ne*ne])
    rotation = _rotation_from_axis(optic_axis_theta_deg, optic_axis_phi_deg)
    epsilon = rotation @ principal @ rotation.T
    kt, kp = np.deg2rad([propagation_theta_deg, propagation_phi_deg])
    direction = np.array([np.sin(kt)*np.cos(kp), np.sin(kt)*np.sin(kp), np.cos(kt)])
    cross = np.array([[0., -direction[2], direction[1]], [direction[2], 0., -direction[0]],
                      [-direction[1], direction[0], 0.]])
    # epsilon E = n^2[-s x mu^-1(s x E)] for exp(i k.r-i omega t).
    eigenvalues, vectors = eig(epsilon, -(cross @ cross)/mu_relative)
    modes = []
    for value, vector in zip(eigenvalues, vectors.T):
        if not np.isfinite(value) or abs(value) > 1e12:
            continue
        index = np.lib.scimath.sqrt(value)
        if index.real < 0 or (abs(index.real) < 1e-14 and index.imag < 0): index = -index
        vector = vector/np.sqrt(np.vdot(vector, vector).real)
        modes.append({"n_complex": {"real": float(index.real), "imag": float(index.imag)},
            "polarization": [{"real": float(x.real), "imag": float(x.imag)} for x in vector],
            "longitudinal_electric_fraction": float(abs(np.dot(direction, vector))**2)})
    modes.sort(key=lambda row: (row["n_complex"]["real"], row["n_complex"]["imag"]))
    if len(modes) != 2:
        raise ValueError("The homogeneous Maxwell eigenproblem did not return two finite transverse modes")
    tensor = [[{"real": float(v.real), "imag": float(v.imag)} for v in row] for row in epsilon]
    return {"epsilon_relative": tensor, "mu_relative": mu_relative, "modes": modes,
        "ordinary_index": {"real": float(no.real), "imag": float(no.imag)},
        "extraordinary_index": {"real": float(ne.real), "imag": float(ne.imag)},
        "contributions": {"thermal_delta_n": float(dn_dT_per_K*(temperature_K-reference_temperature_K)),
            "kerr_delta_n": float(n2_m2_W*intensity_W_m2), "gain_delta_k": float(gain_k)},
        "gain_convention": "positive gain_per_m gives I(z)=I(0) exp(g z) and a negative imaginary refractive index",
        "scope": "Exact local constitutive tensor and homogeneous bulk eigenwaves. This result is not yet a layered FMM, waveguide, thermal-diffusion, gain-saturation, or nonlinear self-consistent solve."}
