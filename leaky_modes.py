"""Complex-frequency poles of normal-incidence uniform multilayer stacks."""

from __future__ import annotations

import numpy as np
from scipy.optimize import root

from materials import material_n
from stack import StackModel


def _indices(model: StackModel, wavelength_um: float, remove_loss: bool) -> list[complex]:
    values = []
    for layer in model.layers:
        if layer.kind != "uniform":
            raise ValueError("The current complex-pole solver supports uniform layers only")
        value = material_n(layer.background_material, wavelength_um, layer.background_n)
        values.append(complex(value.real if remove_loss else value))
    return values


def _denominator(k0: complex, indices: list[complex], thicknesses: list[float],
                 incident_n: float, exit_n: float) -> complex:
    matrix = np.eye(2, dtype=complex)
    for index, thickness in zip(indices, thicknesses):
        phase = k0*index*thickness
        layer = np.asarray([[np.cos(phase), 1j*np.sin(phase)/index],
                            [1j*index*np.sin(phase), np.cos(phase)]])
        matrix = matrix @ layer
    return (incident_n*matrix[0, 0] + incident_n*exit_n*matrix[0, 1]
            + matrix[1, 0] + exit_n*matrix[1, 1])


def _search(model, center_um, initial_q, remove_loss):
    indices = _indices(model, center_um, remove_loss)
    thicknesses = [layer.thickness_um for layer in model.layers]
    initial = 2*np.pi/center_um * (1+1j/(2*initial_q))
    def equations(value):
        denominator = _denominator(complex(value[0], value[1]), indices, thicknesses,
                                   model.incident_n, model.exit_n)
        return [denominator.real, denominator.imag]
    solution = root(equations, [initial.real, initial.imag], method="hybr")
    k0 = complex(*solution.x)
    residual = abs(_denominator(k0, indices, thicknesses, model.incident_n, model.exit_n))
    if not solution.success or residual > 1e-7 or abs(k0.imag) < 1e-14:
        raise ValueError("No stable pole was found near the supplied wavelength and Q guess")
    wavelength = 2*np.pi/k0
    q = abs(k0.real/(2*k0.imag))
    return {"k0_real_rad_um": float(k0.real), "k0_imag_rad_um": float(k0.imag),
            "complex_wavelength_um": {"real": float(wavelength.real),
                                      "imag": float(wavelength.imag)},
            "resonance_wavelength_um": float(2*np.pi/k0.real), "Q": float(q),
            "denominator_residual": float(residual), "iterations": int(solution.nfev)}


def solve_leaky_mode(model: StackModel, center_um: float, initial_q: float = 100.) -> dict:
    """Find a frozen-index outgoing pole and estimate radiative/absorptive Q."""
    model.validate()
    if abs(model.theta_deg) > 1e-12 or model.polarization == "unpolarized":
        raise ValueError("The uniform-stack pole solver currently requires normal-incidence coherent s or p input")
    if center_um <= 0 or initial_q <= .5:
        raise ValueError("Use a positive center wavelength and initial Q greater than 0.5")
    total = _search(model, center_um, initial_q, False)
    lossless = _search(model, center_um, initial_q, True)
    inverse_abs = max(0., 1/total["Q"]-1/lossless["Q"])
    absorptive_q = None if inverse_abs <= 1e-14 else float(1/inverse_abs)
    return {"total": total, "lossless_radiative": lossless,
            "absorptive_Q_from_inverse_Q_difference": absorptive_q,
            "index_model": "Each layer index is evaluated at the real search wavelength and frozen during complex continuation.",
            "scope": "Normal-incidence poles of isotropic uniform multilayer stacks. Patterned FMM poles, dispersive analytic continuation, and quasinormal-mode fields are not included.",
            "convention": "The transfer-matrix time convention may place decaying poles at positive Im(k0); Q uses |Re(k0)/(2 Im(k0))|.",
            "checks": ["Repeat from several center wavelengths and Q guesses to identify the same pole.",
                       "Require a small denominator residual and convergence of the pole when layers or material data change.",
                       "Use the lossless rerun as a radiative-Q estimate only within the frozen-index model."]}
