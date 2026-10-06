"""Transparent one and two parameter sweeps for the finite-stack solver."""

from __future__ import annotations

from dataclasses import asdict, replace

import numpy as np
import run_jobs

from stack import Layer, StackModel, solve_stack, stack_convergence


GLOBAL_PARAMETERS = {
    "wavelength_um", "theta_deg", "phi_deg", "period_x_um", "period_y_um",
    "lattice_angle_deg", "incident_n", "exit_n", "order_budget", "grid_size",
}
LAYER_PARAMETERS = {"thickness_um", "background_n", "feature_n", "fill_x",
                    "fill_y", "inner_radius", "offset_x", "offset_y"}


def parameter_label(path: str) -> str:
    if path in GLOBAL_PARAMETERS:
        return path
    parts = path.split(".")
    if len(parts) == 3 and parts[0] == "layer" and parts[1].isdigit() \
            and parts[2] in LAYER_PARAMETERS:
        return f"layer {int(parts[1])+1} {parts[2]}"
    raise ValueError("Unsupported sweep parameter")


def set_parameter(model: StackModel, path: str, value: float) -> StackModel:
    if not np.isfinite(value):
        raise ValueError("Sweep values must be finite")
    parameter_label(path)
    if path in GLOBAL_PARAMETERS:
        if path in ("order_budget", "grid_size"):
            value = int(round(value))
        candidate = replace(model, **{path: value})
    else:
        _, raw_index, name = path.split(".")
        index = int(raw_index)
        if not 0 <= index < len(model.layers):
            raise ValueError("Sweep layer number is outside the current stack")
        layers = list(model.layers)
        layers[index] = replace(layers[index], **{name: value})
        candidate = replace(model, layers=tuple(layers))
    candidate.validate()
    return candidate


def sweep_values(start: float, stop: float, points: int) -> np.ndarray:
    if not np.isfinite([start, stop]).all() or not 2 <= points <= 201 or start >= stop:
        raise ValueError("Use 2–201 points with finite increasing sweep limits")
    return np.linspace(start, stop, points)


def evaluate_point(model: StackModel, path_x: str, value_x: float,
                   quantity: str, path_y: str | None = None,
                   value_y: float | None = None) -> dict:
    if quantity not in ("R", "T", "A", "R0", "T0"):
        raise ValueError("Sweep quantity must be R, T, A, R0, or T0")
    current = set_parameter(model, path_x, value_x)
    if path_y:
        if path_y == path_x:
            raise ValueError("Choose two different sweep parameters")
        current = set_parameter(current, path_y, float(value_y))
    result = solve_stack(current)
    return {"x": float(value_x), "y": None if not path_y else float(value_y),
            "value": float(result[quantity]), "quantity": quantity,
            "R": result["R"], "T": result["T"], "A": result["A"],
            "R0": result["R0"], "T0": result["T0"],
            "actual_orders": result["actual_orders"]}


def run_sweep(model: StackModel, path_x: str, start_x: float, stop_x: float,
              points_x: int, quantity: str, path_y: str | None = None,
              start_y: float | None = None, stop_y: float | None = None,
              points_y: int | None = None) -> dict:
    """Synchronous reference implementation used by tests and scripted users."""
    xs = sweep_values(start_x, stop_x, points_x)
    ys = [None] if not path_y else sweep_values(float(start_y), float(stop_y), int(points_y))
    if len(xs)*len(ys) > 40401:
        raise ValueError("A sweep may contain at most 40401 points")
    rows = []
    for y in ys:
        for x in xs:
            try:
                rows.append({**evaluate_point(model, path_x, x, quantity, path_y, y),
                             "status": "solved"})
            except (ValueError, np.linalg.LinAlgError) as exc:
                rows.append({"x": float(x), "y": None if y is None else float(y),
                             "value": None, "status": str(exc)})
    solved = [row for row in rows if row["value"] is not None]
    if not solved:
        raise ValueError("No sweep points could be solved")
    best_max = max(solved, key=lambda row: row["value"])
    best_min = min(solved, key=lambda row: row["value"])
    return {"x_parameter": path_x, "x_label": parameter_label(path_x),
            "x_values": xs.tolist(), "y_parameter": path_y,
            "y_label": None if not path_y else parameter_label(path_y),
            "y_values": [] if not path_y else np.asarray(ys).tolist(),
            "quantity": quantity, "rows": rows,
            "best_maximum": best_max, "best_minimum": best_min,
            "model": asdict(model)}


def diffraction_order_sweep(model: StackModel, path: str, start: float,
                            stop: float, points: int) -> dict:
    """Track every propagating diffraction order while one parameter changes.

    An order missing from a successfully solved point is evanescent (zero
    far-field power) at that point. Failed points remain explicit gaps so a
    Rayleigh cutoff or numerical failure is never drawn as a physical zero.
    """
    if model.polarization == "unpolarized":
        raise ValueError("Choose coherent s or p polarization to track diffraction-order amplitudes and phases")
    xs = sweep_values(start, stop, points)
    rows = []
    order_keys: set[tuple[int, int]] = set()
    for index, value in enumerate(xs):
        run_jobs.progress(index, len(xs), "Separating diffraction orders")
        try:
            current = set_parameter(model, path, float(value))
            result = solve_stack(current)
            orders = [{key: item[key] for key in
                       ("m", "n", "R", "T", "reflected_angle_deg", "transmitted_angle_deg")}
                      for item in result.get("order_amplitudes", [])]
            order_keys.update((item["m"], item["n"]) for item in orders)
            rows.append({"x": float(value), "status": "solved", "R": result["R"],
                         "T": result["T"], "A": result["A"],
                         "actual_orders": result["actual_orders"], "orders": orders})
        except (ValueError, np.linalg.LinAlgError) as exc:
            rows.append({"x": float(value), "status": str(exc), "orders": []})
    run_jobs.progress(len(xs), len(xs), "Validating the midpoint")
    solved = [row for row in rows if row["status"] == "solved"]
    if not solved:
        causes = sorted({row["status"] for row in rows})
        raise ValueError("No diffraction-sweep points could be solved. " + " | ".join(causes[:3]))

    center_value = float(xs[len(xs)//2])
    validation = None
    try:
        check = stack_convergence(set_parameter(model, path, center_value))
        validation = {"x": center_value, "converged": check["converged"],
                      "physical_balance_ok": check["physical_balance_ok"],
                      "max_change": check["max_change"], "tolerance": check["tolerance"],
                      "actual_orders": [sample["actual_orders"] for sample in check["samples"]]}
    except (ValueError, np.linalg.LinAlgError) as exc:
        validation = {"x": center_value, "converged": False,
                      "physical_balance_ok": False, "error": str(exc)}
    return {"parameter": path, "parameter_label": parameter_label(path),
            "x_values": xs.tolist(), "rows": rows,
            "orders": [{"m": m, "n": n} for m, n in sorted(order_keys)],
            "solved_points": len(solved), "requested_points": len(rows),
            "validation": validation, "model": asdict(model),
            "interpretation": "Order efficiencies are fractions of incident power. Only propagating far-field orders carry R or T; an order becomes zero when it is evanescent. Failed points are gaps, not zeros."}
