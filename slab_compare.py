"""Planar slab-mode and grating phase-matching comparison for GMR design."""

from dataclasses import asdict, replace
import math

import numpy as np

from stack import StackModel, _grid, solve_stack
from waveguide import WaveguideModel, solve_waveguide


def slab_phase_match(stack: StackModel, slab: WaveguideModel,
                     maximum_order: int = 2) -> dict:
    """Compare bound slab beta with |k_parallel + m b1 + n b2|."""
    if not 0 <= maximum_order <= 5:
        raise ValueError("Use reciprocal orders from 0 to 5")
    stack.validate()
    solved = solve_waveguide(slab)
    gamma = math.radians(stack.lattice_angle_deg)
    a1 = np.array([stack.period_x_um, 0.0])
    a2 = np.array([stack.period_y_um*math.cos(gamma),
                   stack.period_y_um*math.sin(gamma)])
    direct = np.column_stack((a1, a2))
    reciprocal = 2*np.pi*np.linalg.inv(direct).T
    k0 = 2*np.pi/slab.wavelength_um
    theta, phi = math.radians(stack.theta_deg), math.radians(stack.phi_deg)
    k_parallel = stack.incident_n*k0*math.sin(theta)*np.array([math.cos(phi),
                                                               math.sin(phi)])
    rows = []
    for mode in solved["modes"]:
        beta = mode["beta_rad_um"]
        for m in range(-maximum_order, maximum_order+1):
            for n in range(-maximum_order, maximum_order+1):
                if m == 0 and n == 0:
                    continue
                coupled = k_parallel + m*reciprocal[:, 0] + n*reciprocal[:, 1]
                magnitude = float(np.linalg.norm(coupled))
                mismatch = magnitude-beta
                rows.append({"mode": mode["name"], "mode_order": mode["order"],
                             "m": m, "n": n, "beta_rad_um": beta,
                             "coupling_wavevector_rad_um": magnitude,
                             "mismatch_rad_um": mismatch,
                             "relative_mismatch": mismatch/beta,
                             "absolute_relative_mismatch": abs(mismatch/beta),
                             "n_eff": mode["n_eff"],
                             "confinement_fraction": mode["confinement_fraction"]})
    rows.sort(key=lambda row: row["absolute_relative_mismatch"])
    best_by_mode = []
    for mode in solved["modes"]:
        candidate = next((row for row in rows if row["mode_order"] == mode["order"]), None)
        if candidate:
            best_by_mode.append(candidate)
    return {"slab": solved, "stack": asdict(stack), "maximum_order": maximum_order,
            "incident_k_parallel_rad_um": k_parallel.tolist(),
            "reciprocal_vectors_rad_um": {"b1": reciprocal[:, 0].tolist(),
                                            "b2": reciprocal[:, 1].tolist()},
            "best_by_mode": best_by_mode, "rows": rows,
            "equation": "beta ≈ |k_parallel + m b1 + n b2|",
            "warning": "Momentum matching is a design estimate. The periodic perturbation shifts the mode and sets coupling strength, linewidth, and radiation loss; confirm the actual feature with FMM and field/convergence checks."}


def compare_stack_layer(stack: StackModel, layer_index: int, polarization: str,
                        top_n: float, bottom_n: float,
                        maximum_order: int = 2) -> dict:
    """Create an explicitly homogenized three-medium slab from one stack layer."""
    if not 0 <= layer_index < len(stack.layers):
        raise ValueError("Core layer number is outside the finite stack")
    if polarization not in ("TE", "TM"):
        raise ValueError("Slab comparison polarization must be TE or TM")
    if not np.isfinite([top_n, bottom_n]).all() or min(top_n, bottom_n) <= 0:
        raise ValueError("Cladding indices must be finite and positive")
    layer = stack.layers[layer_index]
    mean_epsilon = complex(np.mean(_grid(layer, stack)))
    effective_index = np.sqrt(mean_epsilon)
    if abs(effective_index.imag) > 1e-9:
        raise ValueError("The selected layer has loss; the bound slab comparison is lossless only")
    slab = WaveguideModel(wavelength_um=stack.wavelength_um,
                          thickness_um=layer.thickness_um,
                          core_material="dielectric", core_n=float(effective_index.real),
                          top_material="dielectric", top_n=float(top_n),
                          bottom_material="dielectric", bottom_n=float(bottom_n),
                          polarization=polarization)
    result = slab_phase_match(stack, slab, maximum_order)
    result["homogenization"] = {"layer": layer_index+1,
        "method": "unit-cell arithmetic mean of epsilon",
        "mean_epsilon": float(mean_epsilon.real),
        "effective_core_n": float(effective_index.real),
        "exact_for_uniform_layer": layer.kind == "uniform"}
    return result


def dispersion_comparison(stack: StackModel, layer_index: int, polarization: str,
                          top_n: float, bottom_n: float, maximum_order: int,
                          wavelength_start: float, wavelength_stop: float,
                          points: int, mode_order: int = 0) -> dict:
    """Track a homogenized slab branch beside the actual FMM spectrum."""
    if not np.isfinite([wavelength_start, wavelength_stop]).all() or not 0 < wavelength_start < wavelength_stop:
        raise ValueError("Use positive increasing wavelength limits")
    if not 5 <= points <= 41:
        raise ValueError("Use 5–41 dispersion wavelengths")
    if not 0 <= mode_order <= 49:
        raise ValueError("Mode order must be between 0 and 49")
    rows, failures = [], []
    for wavelength in np.linspace(wavelength_start, wavelength_stop, points):
        current = replace(stack, wavelength_um=float(wavelength))
        try:
            comparison = compare_stack_layer(current, layer_index, polarization,
                                             top_n, bottom_n, maximum_order)
            mode = next((row for row in comparison["best_by_mode"]
                         if row["mode_order"] == mode_order), None)
            scattering = solve_stack(current)
            rows.append({"wavelength_um": float(wavelength),
                         "R": scattering["R"], "T": scattering["T"],
                         "A": scattering["A"],
                         "n_eff": None if mode is None else mode["n_eff"],
                         "relative_mismatch": None if mode is None else mode["relative_mismatch"],
                         "m": None if mode is None else mode["m"],
                         "n": None if mode is None else mode["n"]})
        except (ValueError, np.linalg.LinAlgError) as exc:
            failures.append({"wavelength_um": float(wavelength), "reason": str(exc)})
    valid = [row for row in rows if row["relative_mismatch"] is not None]
    closest = min(valid, key=lambda row: abs(row["relative_mismatch"])) if valid else None
    peak = max(rows, key=lambda row: row["R"]) if rows else None
    return {"rows": rows, "failures": failures, "mode_order": mode_order,
            "polarization": polarization, "closest_phase_match": closest,
            "sampled_reflectance_peak": peak,
            "warning": "The slab branch uses wavelength-dependent material data but a homogenized selected layer. FMM peaks can shift from zero mismatch and may have other physical origins."}
