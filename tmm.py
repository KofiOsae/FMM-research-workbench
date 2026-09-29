"""Scattering-recursion transfer matrix calculation for uniform isotropic stacks."""

from dataclasses import asdict

import numpy as np

from materials import material_n
from stack import StackModel


def _one_polarization(model: StackModel, polarization: str) -> dict:
    n = [complex(model.incident_n)]
    n.extend(material_n(layer.background_material, model.wavelength_um,
                        layer.background_n) for layer in model.layers)
    n.append(complex(model.exit_n))
    kx = model.incident_n * np.sin(np.deg2rad(model.theta_deg))
    z = [np.lib.scimath.sqrt(value * value - kx * kx) for value in n]
    z = [(-value if value.imag < -1e-12 or
          (abs(value.imag) < 1e-12 and value.real < 0) else value) for value in z]
    q = [value if polarization == "s" else index * index / value
         for index, value in zip(n, z)]
    reflection = 0j
    transmission = 1+0j
    for j in range(len(n)-2, -1, -1):
        r_interface = (q[j] - q[j+1]) / (q[j] + q[j+1])
        t_interface = 2*q[j] / (q[j] + q[j+1])
        phase = (np.exp(2j*np.pi*z[j+1]*model.layers[j].thickness_um /
                        model.wavelength_um) if j < len(model.layers) else 1+0j)
        denominator = 1 + r_interface*reflection*phase*phase
        transmission = t_interface*transmission*phase/denominator
        reflection = (r_interface + reflection*phase*phase)/denominator
    R = float(abs(reflection)**2)
    T = float(np.real(q[-1])/np.real(q[0]) * abs(transmission)**2)
    A = float(1-R-T)
    return {"R": R, "T": T, "A": 0.0 if abs(A) < 1e-12 else A,
            "r_phase_deg": float(np.angle(reflection, deg=True)) if abs(reflection)>1e-10 else None,
            "t_phase_deg": float(np.angle(transmission, deg=True)) if abs(transmission)>1e-10 else None,
            "r": {"real":float(reflection.real), "imag":float(reflection.imag)},
            "t": {"real":float(transmission.real), "imag":float(transmission.imag)}}


def solve_tmm(model: StackModel) -> dict:
    """Return power fractions for a laterally uniform stack; no diffraction orders."""
    model.validate()
    if any(layer.kind != "uniform" for layer in model.layers):
        raise ValueError("Transfer matrix requires every layer to be uniform. Use FMM for patterned layers.")
    if model.polarization == "unpolarized":
        s, p = (_one_polarization(model, key) for key in ("s", "p"))
        result = {key: (s[key]+p[key])/2 for key in ("R", "T", "A")}
        result.update(r_phase_deg=None, t_phase_deg=None, polarization_channels={"s":s,"p":p})
    else:
        result = _one_polarization(model, model.polarization)
    return {**result, "method": "uniform-stack scattering recursion",
            "phase_convention": "Wrapped degrees of tangential electric-field amplitudes; reference planes are the first and last interfaces, exp(-i omega t). Null at vanishing amplitude or for incoherent unpolarized input; individual s/p channels remain available.",
            "model": asdict(model)}
