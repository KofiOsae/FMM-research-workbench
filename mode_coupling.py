"""Waveguide-port overlap for analytic and imported transverse illumination.

The calculation takes place on one cross-sectional port plane.  It does not
propagate through a finite grating; its result is the fraction of a supplied
complex port field that matches each solved bound waveguide mode.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from vector_modes import solve_vector_modes


def _complex_field(fields: dict, name: str) -> np.ndarray:
    return (np.asarray(fields[name + "_real"], dtype=float)
            + 1j*np.asarray(fields[name + "_imag"], dtype=float))


def _integral(values: np.ndarray, dx: float, dy: float) -> complex:
    return np.sum(values)*dx*dy


def _gaussian_source(x: np.ndarray, y: np.ndarray, wavelength_um: float,
                     settings: dict) -> dict[str, np.ndarray]:
    wx = float(settings.get("waist_x_um", .8))
    wy = float(settings.get("waist_y_um", .8))
    x0 = float(settings.get("center_x_um", 0))
    y0 = float(settings.get("center_y_um", 0))
    rotation = np.deg2rad(float(settings.get("rotation_deg", 0)))
    order = float(settings.get("super_gaussian_order", 1))
    angle = np.deg2rad(float(settings.get("polarization_angle_deg", 0)))
    relative_phase = np.deg2rad(float(settings.get("polarization_phase_deg", 0)))
    theta = np.deg2rad(float(settings.get("theta_deg", 0)))
    phi = np.deg2rad(float(settings.get("phi_deg", 0)))
    medium_n = float(settings.get("medium_n", 1))
    values = np.asarray([wx, wy, x0, y0, rotation, order, angle,
                         relative_phase, theta, phi, medium_n], dtype=float)
    if not np.isfinite(values).all() or wx <= 0 or wy <= 0 or order < 1 or medium_n <= 0:
        raise ValueError("Beam waists and medium index must be positive; super-Gaussian order must be at least 1")
    if abs(np.rad2deg(theta)) >= 80:
        raise ValueError("Use a port-plane beam tilt below 80 degrees")
    xx, yy = np.meshgrid(x, y)
    xr = np.cos(rotation)*(xx-x0)+np.sin(rotation)*(yy-y0)
    yr = -np.sin(rotation)*(xx-x0)+np.cos(rotation)*(yy-y0)
    envelope = np.exp(-((xr/wx)**2+(yr/wy)**2)**order)
    kx = medium_n*np.sin(theta)*np.cos(phi)
    ky = medium_n*np.sin(theta)*np.sin(phi)
    kz = medium_n*np.cos(theta)
    phase = np.exp(1j*2*np.pi/wavelength_um*(kx*xx+ky*yy))
    ex = np.cos(angle)*envelope*phase
    ey = np.sin(angle)*np.exp(1j*relative_phase)*envelope*phase
    ez = -(kx*ex+ky*ey)/kz
    # Relative Maxwell units used by the vector modesolver: H = n k-hat × E.
    khat = np.asarray([kx, ky, kz])/medium_n
    electric = np.stack((ex, ey, ez), axis=-1)
    magnetic = medium_n*np.cross(np.broadcast_to(khat, electric.shape), electric)
    return {"Ex": ex, "Ey": ey, "Ez": ez,
            "Hx": magnetic[..., 0], "Hy": magnetic[..., 1],
            "Hz": magnetic[..., 2]}


def _imported_source(x_target: np.ndarray, y_target: np.ndarray,
                     imported: dict) -> tuple[dict[str, np.ndarray], list[str]]:
    x = np.asarray(imported.get("x_um", []), dtype=float)
    y = np.asarray(imported.get("y_um", []), dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) < 2 or len(y) < 2:
        raise ValueError("Imported fields need increasing x_um and y_um axes with at least two values each")
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.all(np.diff(x) > 0) or not np.all(np.diff(y) > 0):
        raise ValueError("Imported x_um and y_um axes must be finite and strictly increasing")
    xx, yy = np.meshgrid(x_target, y_target)
    points = np.column_stack((yy.ravel(), xx.ravel()))
    result, present = {}, []
    for name in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"):
        real_key, imag_key = name + "_real", name + "_imag"
        if real_key not in imported:
            continue
        real = np.asarray(imported[real_key], dtype=float)
        imag = np.asarray(imported.get(imag_key, np.zeros_like(real)), dtype=float)
        if real.shape != (len(y), len(x)) or imag.shape != real.shape or not np.isfinite(real).all() or not np.isfinite(imag).all():
            raise ValueError(f"Imported {name} arrays must have shape len(y_um) × len(x_um) and finite values")
        field = real+1j*imag
        interpolation = RegularGridInterpolator((y, x), field, bounds_error=False,
                                                  fill_value=0)
        result[name] = interpolation(points).reshape(len(y_target), len(x_target))
        present.append(name)
    if "Ex" not in result or "Ey" not in result:
        raise ValueError("Imported illumination must include Ex_real and Ey_real; imaginary columns are optional")
    for name in ("Ez", "Hx", "Hy", "Hz"):
        result.setdefault(name, np.zeros_like(result["Ex"]))
    return result, present


def _electric_overlap(source: dict[str, np.ndarray], mode: dict[str, np.ndarray],
                      dx: float, dy: float) -> tuple[complex, float]:
    dot = sum(source[name]*np.conj(mode[name]) for name in ("Ex", "Ey", "Ez"))
    source_norm = float(np.real(_integral(sum(abs(source[name])**2 for name in ("Ex", "Ey", "Ez")), dx, dy)))
    mode_norm = float(np.real(_integral(sum(abs(mode[name])**2 for name in ("Ex", "Ey", "Ez")), dx, dy)))
    if source_norm <= 1e-30 or mode_norm <= 1e-30:
        raise ValueError("Source or mode has zero electric-field norm on the port")
    amplitude = _integral(dot, dx, dy)/np.sqrt(source_norm*mode_norm)
    return amplitude, float(abs(amplitude)**2)


def _power_overlaps(source: dict[str, np.ndarray], mode: dict[str, np.ndarray],
                    dx: float, dy: float) -> dict | None:
    source_power = .5*float(np.real(_integral(
        source["Ex"]*np.conj(source["Hy"])-source["Ey"]*np.conj(source["Hx"]), dx, dy)))
    mode_power = .5*float(np.real(_integral(
        mode["Ex"]*np.conj(mode["Hy"])-mode["Ey"]*np.conj(mode["Hx"]), dx, dy)))
    if source_power <= 1e-20 or mode_power <= 1e-20:
        return None
    first = _integral(source["Ex"]*np.conj(mode["Hy"])
                      -source["Ey"]*np.conj(mode["Hx"]), dx, dy)
    second = _integral(np.conj(mode["Ex"])*source["Hy"]
                       -np.conj(mode["Ey"])*source["Hx"], dx, dy)
    scale = 4*np.sqrt(source_power*mode_power)
    forward = (first+second)/scale
    backward = (first-second)/scale
    return {"source_power_relative": source_power,
            "mode_power_relative": mode_power,
            "forward_amplitude": {"real": float(forward.real), "imag": float(forward.imag)},
            "backward_amplitude": {"real": float(backward.real), "imag": float(backward.imag)},
            "forward_efficiency": float(abs(forward)**2),
            "backward_efficiency": float(abs(backward)**2)}


def _project_on_solved(solved: dict, mode_solver: dict, source: dict,
                       include_maps: bool = True) -> dict:
    """Project a source on one already-solved, uniformly normalized port basis."""
    x, y = np.asarray(solved["x_um"]), np.asarray(solved["y_um"])
    source_type = str(source.get("type", "gaussian"))
    if source_type in ("gaussian", "super_gaussian"):
        settings = dict(source)
        if source_type == "gaussian":
            settings["super_gaussian_order"] = 1
        source_fields = _gaussian_source(x, y, float(mode_solver["wavelength_um"]), settings)
        present = ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]
        interpolation = None
    elif source_type == "imported":
        source_fields, present = _imported_source(x, y, source.get("field", {}))
        interpolation = "linear interpolation to the waveguide cell centers; values outside the imported rectangle are zero"
    else:
        raise ValueError("Illumination type must be gaussian, super_gaussian, or imported")
    has_power_fields = "Hx" in present and "Hy" in present
    dx, dy = solved["mesh"]["dx_um"], solved["mesh"]["dy_um"]
    rows, mode_fields = [], []
    for item in solved["modes"]:
        mode = {name: _complex_field(item["fields"], name)
                for name in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")}
        mode_fields.append(mode)
        amplitude, electric_efficiency = _electric_overlap(source_fields, mode, dx, dy)
        power = _power_overlaps(source_fields, mode, dx, dy) if has_power_fields else None
        rows.append({"mode": item["mode"], "n_eff": item["n_eff"],
                     "te_like_fraction": item["te_like_fraction"],
                     "core_electric_fraction": item["core_electric_fraction"],
                     "electric_overlap_amplitude": {"real": float(amplitude.real),
                                                      "imag": float(amplitude.imag)},
                     "electric_profile_overlap": electric_efficiency,
                     "power_overlap": power})
    best_key = ("forward_efficiency" if has_power_fields else None)
    best = max(rows, key=lambda row: row["power_overlap"][best_key]
               if best_key and row["power_overlap"] else row["electric_profile_overlap"])
    intensity = sum(abs(source_fields[name])**2 for name in ("Ex", "Ey", "Ez"))
    maximum = max(float(np.max(intensity)), 1e-300)
    result = {"source_type": source_type, "source_components": present,
            "x_um": x.tolist(), "y_um": y.tolist(),
            "modes": rows, "best_mode": best["mode"],
            "mesh": solved["mesh"], "geometry": solved["geometry"],
            "interpolation": interpolation,
            "formalism": {
                "electric_profile": "|integral E_source dot E_mode* dA|^2 / (integral |E_source|^2 dA integral |E_mode|^2 dA)",
                "directional_power": "reciprocity projection from the transverse E and H fields on one uniform port plane",
                "reference": "A. W. Snyder and J. D. Love, Optical Waveguide Theory, mode orthogonality and excitation chapters; vector modes follow Fallahkhair, Li, and Murphy, JLT 26, 1423 (2008), DOI 10.1109/JLT.2008.923643"},
            "quantity_status": {
                "power_overlap_available": has_power_fields,
                "electric_profile_overlap": "Normalized complex electric-field profile and polarization overlap; it is not a propagated device efficiency.",
                "power_overlap": "Reciprocity projection on this port plane. It is a modal coupling efficiency only when the imported/analytic field is a physically normalized field at this same uniform waveguide port."},
            "device_limit": "No finite grating, fiber-to-chip propagation, apodization, substrate leakage, or transition loss is included. Multiply neither overlap by a diffraction efficiency without a validated common port normalization and finite-device model."}
    if has_power_fields:
        forward = float(sum((row["power_overlap"] or {}).get("forward_efficiency", 0)
                            for row in rows))
        backward = float(sum((row["power_overlap"] or {}).get("backward_efficiency", 0)
                             for row in rows))
        result["solved_mode_power_accounting"] = {
            "forward_fraction_in_solved_modes": forward,
            "forward_unresolved_or_radiation_fraction": max(0.0, 1-forward),
            "backward_projection_sum": backward,
            "forward_sum_within_unit_interval": forward <= 1.0+5e-3,
            "note": "The unresolved fraction includes radiation, continuum fields, modes not requested, and numerical error. Backward projection is a direction diagnostic and is not added to the forward budget."}
    if include_maps:
        best_index = next(i for i, row in enumerate(rows) if row["mode"] == best["mode"])
        best_fields = mode_fields[best_index]
        mode_intensity = sum(abs(best_fields[name])**2 for name in ("Ex", "Ey", "Ez"))
        overlap_density = abs(sum(source_fields[name]*np.conj(best_fields[name])
                                  for name in ("Ex", "Ey", "Ez")))
        result["source_E2_normalized"] = (intensity/maximum).real.tolist()
        result["best_mode_E2_normalized"] = (mode_intensity/max(float(np.max(mode_intensity)), 1e-300)).real.tolist()
        result["overlap_density_normalized"] = (overlap_density/max(float(np.max(overlap_density)), 1e-300)).real.tolist()
        result["map_normalization"] = "Each panel is normalized to its own maximum; use the reported overlap and power values for quantitative comparison."
    return result


def mode_port_coupling(mode_solver: dict, source: dict) -> dict:
    """Solve the cross section and project one supplied port field onto its modes."""
    solved = solve_vector_modes(**mode_solver)
    return _project_on_solved(solved, mode_solver, source, include_maps=True)


def mode_port_coupling_sweep(mode_solver: dict, source: dict, parameter: str,
                             start: float, stop: float, points: int) -> dict:
    """Screen one source parameter while reusing one solved waveguide basis."""
    allowed = {"waist_x_um", "waist_y_um", "center_x_um", "center_y_um",
               "theta_deg", "phi_deg", "polarization_angle_deg"}
    if parameter not in allowed:
        raise ValueError("Unsupported port sweep parameter")
    if source.get("type", "gaussian") not in ("gaussian", "super_gaussian"):
        raise ValueError("Automated port sweeps require an analytic Gaussian or super-Gaussian source")
    if not 3 <= int(points) <= 41 or not np.isfinite([start, stop]).all() or stop <= start:
        raise ValueError("Use 3–41 sweep points and increasing finite limits")
    solved = solve_vector_modes(**mode_solver)
    values = np.linspace(float(start), float(stop), int(points))
    rows = []
    for value in values:
        settings = {**source, parameter: float(value)}
        projected = _project_on_solved(solved, mode_solver, settings, include_maps=False)
        best = next(row for row in projected["modes"] if row["mode"] == projected["best_mode"])
        rows.append({"value": float(value), "best_mode": projected["best_mode"],
                     "electric_profile_overlap": best["electric_profile_overlap"],
                     "forward_efficiency": None if best["power_overlap"] is None else
                        best["power_overlap"]["forward_efficiency"],
                     "all_modes": [{"mode": row["mode"],
                                    "electric_profile_overlap": row["electric_profile_overlap"],
                                    "forward_efficiency": None if row["power_overlap"] is None else
                                        row["power_overlap"]["forward_efficiency"]}
                                   for row in projected["modes"]]})
    metric = "forward_efficiency" if rows[0]["forward_efficiency"] is not None else "electric_profile_overlap"
    best_row = max(rows, key=lambda row: row[metric])
    return {"parameter": parameter, "start": float(start), "stop": float(stop),
            "points": int(points), "metric": metric, "rows": rows,
            "best": best_row, "mode_basis_reused": True,
            "interpretation": "This is a source-to-uniform-port matching sweep. The waveguide eigenbasis is solved once because geometry and wavelength remain fixed; it is not a finite-grating propagation sweep."}


def mode_port_coupling_validation(mode_solver: dict, source: dict,
                                  neff_tolerance: float = 5e-4,
                                  overlap_tolerance: float = 5e-3,
                                  profile: str = "research") -> dict:
    """Repeat the port projection at half mesh and track the nearest physical mode."""
    profiles = {"exploratory": (5e-3, 2e-2),
                "research": (5e-4, 5e-3),
                "publication": (1e-4, 1e-3)}
    profile = str(profile).lower()
    if profile not in (*profiles, "custom"):
        raise ValueError("Mode validation profile must be exploratory, research, publication, or custom")
    if profile != "custom":
        neff_tolerance, overlap_tolerance = profiles[profile]
    if not (0 < neff_tolerance <= .1 and 0 < overlap_tolerance <= .25):
        raise ValueError("Mode validation tolerances must be positive and physically bounded")
    coarse = mode_port_coupling(mode_solver, source)
    refined_inputs = {**mode_solver, "mesh_um": float(mode_solver["mesh_um"])/2}
    refined = mode_port_coupling(refined_inputs, source)
    coarse_row = next(row for row in coarse["modes"] if row["mode"] == coarse["best_mode"])
    target_neff = coarse_row["n_eff"]["real"]
    target_te = coarse_row["te_like_fraction"]
    refined_row = min(refined["modes"], key=lambda row:
                      abs(row["n_eff"]["real"]-target_neff)+.05*abs(row["te_like_fraction"]-target_te))
    metric = "forward_efficiency" if coarse_row["power_overlap"] is not None else "electric_profile_overlap"
    coarse_value = (coarse_row["power_overlap"][metric] if coarse_row["power_overlap"] is not None
                    else coarse_row[metric])
    refined_value = (refined_row["power_overlap"][metric] if refined_row["power_overlap"] is not None
                     else refined_row[metric])
    neff_change = abs(refined_row["n_eff"]["real"]-target_neff)
    overlap_change = abs(refined_value-coarse_value)
    passed = neff_change <= neff_tolerance and overlap_change <= overlap_tolerance
    return {"coarse": coarse, "refined": refined,
            "acceptance_profile": {"name": profile,
                "effective_index_tolerance": neff_tolerance,
                "overlap_absolute_tolerance": overlap_tolerance,
                "scope": "Half-mesh stability of the tracked uniform-port mode and source overlap; padding and mode-count convergence remain separate checks."},
            "tracked_mode": {"coarse_mode": coarse_row["mode"], "refined_mode": refined_row["mode"],
                             "tracking_basis": "nearest effective index with TE-like-fraction tie breaking"},
            "metric": metric, "coarse_value": coarse_value, "refined_value": refined_value,
            "effective_index_change": neff_change, "overlap_change": overlap_change,
            "tolerances": {"effective_index": neff_tolerance, "overlap": overlap_tolerance},
            "passed": bool(passed),
            "scope": "Mesh refinement certifies the uniform-port eigenmode projection only. Domain padding, mode count, imported-field sampling, and any finite-device scattering require separate checks."}
