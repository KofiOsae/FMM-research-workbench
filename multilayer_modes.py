"""Lossless planar multilayer modes from a one-dimensional finite-volume eigenproblem."""

from dataclasses import asdict
import math

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

from stack import StackModel, _grid


def solve_multilayer_modes(model: StackModel, polarization: str = "TE",
                           grid_points: int = 1201, modes: int = 8) -> dict:
    if polarization not in ("TE", "TM"):
        raise ValueError("Multilayer polarization must be TE or TM")
    if not 401 <= grid_points <= 4001 or grid_points % 2 == 0:
        raise ValueError("Use an odd vertical grid from 401 to 4001 points")
    if not 1 <= modes <= 20:
        raise ValueError("Request 1–20 multilayer modes")
    model.validate()
    layer_eps = []
    for layer in model.layers:
        value = complex(np.mean(_grid(layer, model)))
        if abs(value.imag) > 1e-9:
            raise ValueError("The multilayer eigenmode solver requires lossless materials")
        layer_eps.append(float(value.real))
    top_eps, bottom_eps = model.incident_n**2, model.exit_n**2
    thicknesses = np.asarray([layer.thickness_um for layer in model.layers])
    total = float(thicknesses.sum())
    padding = max(3*model.wavelength_um, 2*total)
    z = np.linspace(-padding, total+padding, grid_points)
    dz = float(z[1]-z[0])
    edges = np.concatenate(([0.0], np.cumsum(thicknesses)))
    eps = np.full(grid_points, top_eps)
    eps[z > total] = bottom_eps
    for i, value in enumerate(layer_eps):
        eps[(z >= edges[i]) & (z <= edges[i+1])] = value
    interior_eps = eps[1:-1]
    count = len(interior_eps)
    k0 = 2*np.pi/model.wavelength_um
    if polarization == "TE":
        main = -2*np.ones(count)/dz**2 + k0**2*interior_eps
        off = np.ones(count-1)/dz**2
        operator = diags((off, main, off), (-1, 0, 1), format="csc")
        mass = None
    else:
        inv_eps = 1/eps
        interface = 2*inv_eps[:-1]*inv_eps[1:]/(inv_eps[:-1]+inv_eps[1:])
        left, right = interface[:-1], interface[1:]
        main = -(left+right)/dz**2 + k0**2*np.ones(count)
        operator = diags((left[1:]/dz**2, main, right[:-1]/dz**2),
                         (-1, 0, 1), format="csc")
        mass = diags(1/interior_eps, 0, format="csc")
    requested = min(max(modes+6, 12), count-2)
    values, vectors = eigsh(operator, k=requested, M=mass, which="LA")
    order = np.argsort(values)[::-1]
    lower = max(model.incident_n, model.exit_n)
    upper = math.sqrt(max(layer_eps))
    found = []
    for value, vector in zip(values[order], vectors[:, order].T):
        if value <= 0:
            continue
        beta = math.sqrt(float(value))
        neff = beta/k0
        if not lower+1e-5 < neff < upper-1e-5:
            continue
        full = np.zeros(grid_points)
        full[1:-1] = vector
        intensity = full**2
        intensity /= max(float(intensity.max()), 1e-30)
        core_mask = (z >= 0) & (z <= total)
        norm_weight = np.ones_like(z) if polarization == "TE" else 1/eps
        confinement = float(np.trapezoid(intensity[core_mask]*norm_weight[core_mask], z[core_mask])
                            / np.trapezoid(intensity*norm_weight, z))
        found.append({"order": len(found), "name": f"{polarization}{len(found)}",
                      "n_eff": neff, "beta_rad_um": beta,
                      "confinement_fraction": confinement,
                      "z_um": z.tolist(), "profile": intensity.tolist(),
                      "profile_quantity": "|E_y|²" if polarization == "TE" else "|H_y|²"})
        if len(found) >= modes:
            break
    return {"model": asdict(model), "polarization": polarization,
            "modes": found, "z_um": z.tolist(), "epsilon": eps.tolist(),
            "homogenized_layer_epsilon": layer_eps, "grid_points": grid_points,
            "padding_um": padding,
            "note": "Every finite layer is homogenized by unit-cell mean epsilon. Dirichlet boundaries are placed in padded claddings; repeat with a finer grid and larger-domain formulation for publication."}
