"""Closed-loop polarization winding from background-subtracted resonance radiation."""

from __future__ import annotations

from dataclasses import asdict
import math
import numpy as np

from stack import StackModel, solve_stack, stack_convergence


def _complex(value: dict) -> complex:
    return complex(float(value["real"]), float(value["imag"]))


def _charge_from_angles(angles_rad: np.ndarray) -> tuple[float, np.ndarray]:
    """Winding of a headless polarization axis using its doubled angle."""
    closed = np.r_[2*np.asarray(angles_rad), 2*angles_rad[0]]
    unwrapped = np.unwrap(closed)
    return float((unwrapped[-1]-unwrapped[0])/(4*np.pi)), unwrapped/2


def _model_at_uv(model: StackModel, wavelength: float, u: float, v: float) -> StackModel:
    rho = math.hypot(u, v)
    if rho >= math.sin(math.radians(89)):
        raise ValueError("Loop extends beyond the supported propagating incident light cone")
    return StackModel(**{**asdict(model), "layers": model.layers,
        "wavelength_um": float(wavelength),
        "theta_deg": math.degrees(math.asin(rho)),
        "phi_deg": 0.0 if rho < 1e-14 else math.degrees(math.atan2(v, u))})


def _cartesian_channel(model: StackModel, wavelength: float, u: float, v: float,
                       m: int, n: int, port: str) -> tuple[complex, complex, float]:
    result = solve_stack(_model_at_uv(model, wavelength, u, v))
    row = next((item for item in result.get("order_amplitudes", [])
                if item["m"] == m and item["n"] == n), None)
    if row is None:
        raise ValueError(f"Order ({m},{n}) is not propagating on the complete loop")
    field = row[port]
    return _complex(field["Ex"]), _complex(field["Ey"]), float(row["R" if port == "reflected" else "T"])


def _loop(model: StackModel, center: float, linewidth: float, radius: float,
          loop_points: int, sideband_span: float, search_span: float,
          quantity: str, extremum: str, m: int, n: int, port: str) -> dict:
    rows, previous_center = [], center
    for index in range(loop_points):
        phase = 2*np.pi*index/loop_points
        u, v = radius*np.cos(phase), radius*np.sin(phase)
        scan = np.linspace(previous_center-search_span*linewidth,
                           previous_center+search_span*linewidth, 7)
        candidates = []
        for wavelength in scan:
            current = solve_stack(_model_at_uv(model, wavelength, u, v))
            candidates.append((wavelength, float(current[quantity])))
        values = np.asarray([item[1] for item in candidates])
        pick = int(np.argmax(values) if extremum == "max" else np.argmin(values))
        local_center = float(candidates[pick][0])
        previous_center = local_center
        center_field = _cartesian_channel(model, local_center, u, v, m, n, port)
        estimates = []
        for span in (sideband_span, sideband_span+1):
            low = _cartesian_channel(model, local_center-span*linewidth, u, v, m, n, port)
            high = _cartesian_channel(model, local_center+span*linewidth, u, v, m, n, port)
            background = ((low[0]+high[0])/2, (low[1]+high[1])/2)
            estimates.append(np.asarray([center_field[0]-background[0],
                                         center_field[1]-background[1]]))
        field = estimates[0]
        amplitude = float(np.linalg.norm(field))
        stability = float(np.linalg.norm(estimates[0]-estimates[1])/max(amplitude, 1e-30))
        s1 = abs(field[0])**2-abs(field[1])**2
        s2 = 2*np.real(field[0]*np.conj(field[1]))
        s3 = -2*np.imag(field[0]*np.conj(field[1]))
        s0 = max(float(abs(field[0])**2+abs(field[1])**2), 1e-30)
        orientation = .5*np.arctan2(s2, s1)
        ellipticity = .5*np.arcsin(np.clip(s3/s0, -1, 1))
        convergence = stack_convergence(_model_at_uv(model, local_center, u, v))
        rows.append({"index": index, "loop_angle_deg": float(np.degrees(phase)),
            "u": float(u), "v": float(v), "resonance_wavelength_um": local_center,
            "tracked_quantity": float(candidates[pick][1]),
            "Ex": {"real": float(field[0].real), "imag": float(field[0].imag)},
            "Ey": {"real": float(field[1].real), "imag": float(field[1].imag)},
            "amplitude": amplitude, "orientation_deg": float(np.degrees(orientation)),
            "ellipticity_deg": float(np.degrees(ellipticity)),
            "sideband_relative_change": stability,
            "fourier_converged": bool(convergence["converged"] and convergence["physical_balance_ok"])})
    angles = np.radians([row["orientation_deg"] for row in rows])
    charge, unwrapped = _charge_from_angles(angles)
    for row, angle in zip(rows, unwrapped[:-1]):
        row["unwrapped_orientation_deg"] = float(np.degrees(angle))
    return {"radius": radius, "rows": rows, "charge": charge,
            "minimum_residual_amplitude": min(row["amplitude"] for row in rows),
            "maximum_sideband_relative_change": max(row["sideband_relative_change"] for row in rows),
            "maximum_absolute_ellipticity_deg": max(abs(row["ellipticity_deg"]) for row in rows),
            "all_fourier_converged": all(row["fourier_converged"] for row in rows),
            "maximum_branch_step_um": max(abs(rows[(i+1) % len(rows)]["resonance_wavelength_um"]-rows[i]["resonance_wavelength_um"]) for i in range(len(rows)))}


def polarization_winding(model: StackModel, center_um: float, linewidth_um: float,
        radius: float, loop_points: int = 12, sideband_span: float = 3,
        search_span: float = 2, quantity: str = "R", extremum: str = "max",
        m: int = 0, n: int = 0, port: str = "reflected") -> dict:
    model.validate()
    if model.polarization not in ("s", "p"):
        raise ValueError("Choose coherent s or p incidence")
    if quantity not in ("R", "T", "A") or extremum not in ("max", "min"):
        raise ValueError("Choose R, T, or A and maximum or minimum tracking")
    if port not in ("reflected", "transmitted"):
        raise ValueError("Choose reflected or transmitted radiation")
    if not np.isfinite([center_um, linewidth_um, radius, sideband_span, search_span]).all() or center_um <= 0 or linewidth_um <= 0:
        raise ValueError("Center and linewidth must be positive finite values")
    if not .002 <= radius <= .7 or not 8 <= loop_points <= 32 or sideband_span < 1.5 or not 1 <= search_span <= 5:
        raise ValueError("Use radius 0.002–0.7, 8–32 loop points, sidebands ≥1.5 linewidths, and search span 1–5")
    if center_um-(search_span+sideband_span+1)*linewidth_um <= 0:
        raise ValueError("The search or sideband wavelength reaches zero")
    primary = _loop(model, center_um, linewidth_um, radius, loop_points,
                    sideband_span, search_span, quantity, extremum, m, n, port)
    comparison = _loop(model, center_um, linewidth_um, radius*1.15, loop_points,
                       sideband_span, search_span, quantity, extremum, m, n, port)
    rounded = int(round(primary["charge"]))
    stable_radius = abs(primary["charge"]-comparison["charge"]) < .15
    integer_like = abs(primary["charge"]-rounded) < .15
    linear_enough = primary["maximum_absolute_ellipticity_deg"] <= 10
    sideband_stable = primary["maximum_sideband_relative_change"] <= .25
    branch_stable = primary["maximum_branch_step_um"] <= 2.1*search_span*linewidth_um
    certified = bool(stable_radius and integer_like and linear_enough and sideband_stable and
                     branch_stable and primary["all_fourier_converged"] and
                     primary["minimum_residual_amplitude"] > 1e-8)
    return {"primary": primary, "radius_comparison": comparison,
        "rounded_charge": rounded if certified else None, "certified": certified,
        "checks": {"radius_stable": stable_radius, "integer_like": integer_like,
            "near_linear_polarization": linear_enough, "sideband_stable": sideband_stable,
            "branch_stable": branch_stable,
            "all_fourier_converged": primary["all_fourier_converged"],
            "residual_nonzero": primary["minimum_residual_amplitude"] > 1e-8},
        "basis": "fixed laboratory Cartesian Ex,Ey",
        "method": "closed-loop winding of background-subtracted resonant radiation with wavelength tracking",
        "limitation": "Sideband subtraction approximates an isolated pole residue. A rigorous ideal-BIC proof still requires a patterned outgoing eigenmode or S-matrix pole solver."}
