"""Observable-aware validation, linked sweeps, and diffraction-order maps."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json

import numpy as np

import run_jobs
from observables import observable_spec, observable_value
from stack import StackModel, grid_convergence, prepare, solve_stack, stack_convergence
from sweep import parameter_label, set_parameter, sweep_values


def settings_fingerprint(value: dict) -> str:
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def diffraction_order_map(model: StackModel) -> dict:
    """Classify retained reciprocal orders at both homogeneous ports."""
    obj = prepare(model)
    k0 = 2*np.pi/model.wavelength_um
    rows = []
    for index, (m, n) in enumerate(obj.G):
        kx, ky = float(np.real(obj.kx[index])), float(np.real(obj.ky[index]))
        item = {"m": int(m), "n": int(n), "kx_over_k0": kx/k0,
                "ky_over_k0": ky/k0, "ports": {}}
        for name, refractive_index in (("reflected", model.incident_n),
                                        ("transmitted", model.exit_n)):
            kz2 = (k0*refractive_index)**2-kx*kx-ky*ky
            propagating = bool(kz2 > 1e-12)
            grazing = bool(abs(kz2) <= 1e-12)
            angle = None
            if propagating:
                angle = float(np.degrees(np.arctan2(np.hypot(kx, ky), np.sqrt(kz2))))
            item["ports"][name] = {"status": "grazing" if grazing else
                "propagating" if propagating else "evanescent",
                "angle_deg": angle, "normalized_kz_squared": float(kz2/k0**2)}
        rows.append(item)
    return {"orders": rows, "actual_orders": int(obj.nG),
            "wavelength_um": model.wavelength_um,
            "incident_n": model.incident_n, "exit_n": model.exit_n,
            "interpretation": "Orders inside a port light cone propagate and have an output angle. Orders outside it are evanescent. A grazing order lies at a Rayleigh cutoff and should be shifted slightly before numerical use."}


def adaptive_observable_convergence(model: StackModel, observables: list[dict],
                                    tolerance: float = .01) -> dict:
    if not 1e-8 <= tolerance <= .25:
        raise ValueError("Convergence tolerance must be from 1e-8 to 0.25")
    center = model.order_budget
    budgets = sorted(set((max(3, center-16), center, min(201, center+24))))
    if len(budgets) < 3:
        budgets = [3, 101, 201]
    attempts = []
    while True:
        report = stack_convergence(model, budgets=budgets[-3:], tolerance=tolerance,
                                   observables=observables)
        attempts.append({"budgets": budgets[-3:], "max_change": report["max_change"],
                         "converged": report["converged"]})
        if report["converged"] or budgets[-1] >= 201:
            report["adaptive_attempts"] = attempts
            report["reached_maximum_budget"] = bool(budgets[-1] >= 201)
            return report
        candidate = min(201, budgets[-1]+24)
        if candidate in budgets:
            report["adaptive_attempts"] = attempts
            report["reached_maximum_budget"] = True
            return report
        budgets.append(candidate)


VALIDATION_PROFILES = {
    "exploratory": {"observable_tolerance": 1e-2, "energy_tolerance": 1e-3,
                    "purpose": "Fast screening; insufficient for a publication claim."},
    "research": {"observable_tolerance": 1e-3, "energy_tolerance": 1e-4,
                 "purpose": "Research iteration and candidate selection."},
    "publication": {"observable_tolerance": 2e-4, "energy_tolerance": 2e-5,
                    "purpose": "Strict numerical screen for power and diffraction-order observables."},
}


def validation_report(model: StackModel, observables: list[dict] | None = None,
                      tolerance: float = 1e-3, adaptive: bool = True,
                      profile: str = "research") -> dict:
    model.validate()
    profile = str(profile).lower()
    if profile not in (*VALIDATION_PROFILES, "custom"):
        raise ValueError("Validation profile must be exploratory, research, publication, or custom")
    if profile != "custom":
        tolerance = VALIDATION_PROFILES[profile]["observable_tolerance"]
        energy_tolerance = VALIDATION_PROFILES[profile]["energy_tolerance"]
        purpose = VALIDATION_PROFILES[profile]["purpose"]
    else:
        energy_tolerance = max(1e-10, tolerance/10)
        purpose = "User-defined numerical screen. The exported tolerance is part of the method record."
    requested = [observable_spec(item) | {"tolerance": float(
        item.get("tolerance", tolerance) if profile == "custom" else tolerance)}
                 for item in (observables or [])]
    fourier = (adaptive_observable_convergence(model, requested, tolerance) if adaptive else
               stack_convergence(model, tolerance=tolerance, observables=requested))
    high_budget = int(fourier["samples"][-1]["requested_budget"])
    refined = StackModel(**{**asdict(model), "layers": model.layers,
                            "order_budget": high_budget})
    try:
        grid = grid_convergence(refined, tolerance=tolerance, observables=requested)
        grid_error = None
    except (ValueError, np.linalg.LinAlgError) as exc:
        grid, grid_error = None, str(exc)
    result = solve_stack(refined)
    balance = abs(result["R"]+result["T"]+result["A"]-1)
    material = {"status": "valid", "message":
                "Every selected material supplied optical constants at the submitted wavelength; extrapolation is disabled."}
    checks = {"energy_balance": balance <= energy_tolerance,
              "fourier": bool(fourier["converged"] and fourier["physical_balance_ok"]),
              "geometry_grid": bool(grid and grid["converged"]),
              "material_range": True}
    status = "Validated" if all(checks.values()) else "Not converged"
    submitted = asdict(model)
    return {"status": status, "checks": checks, "energy_balance_residual": balance,
            "acceptance_profile": {"name": profile,
                "observable_absolute_tolerance": tolerance,
                "energy_balance_tolerance": energy_tolerance,
                "purpose": purpose,
                "scope": "Applies to dimensionless power and diffraction-order observables. Resonance center, linewidth, Q, field maps, modes, and statistical yield require their own quantity-specific stability tests."},
            "fourier": fourier, "geometry_grid": grid, "geometry_grid_error": grid_error,
            "material_validity": material, "result": result,
            "order_map": diffraction_order_map(refined),
            "submitted_model": submitted,
            "submitted_model_sha256": settings_fingerprint(submitted),
            "interpretation": "Validation follows total power, every retained propagating order, and each requested observable. Passing is a numerical screen for these refinements, not proof of experimental accuracy or publication readiness by itself."}


def linked_observable_sweep(model: StackModel, driver: str, links: list[dict],
                            start: float, stop: float, points: int,
                            observable: dict, goal: str = "max",
                            validate_each: bool = False,
                            tolerance: float = 1e-3) -> dict:
    if goal not in ("max", "min"):
        raise ValueError("Sweep goal must be max or min")
    spec = observable_spec(observable)
    values = sweep_values(start, stop, points)
    if validate_each and points > 51:
        raise ValueError("Convergence-aware linked sweeps are limited to 51 points")
    parameter_label(driver)
    normalized_links = []
    used = {driver}
    for item in links:
        path = str(item["path"])
        scale, offset = float(item.get("scale", 1)), float(item.get("offset", 0))
        parameter_label(path)
        if path in used or not np.isfinite([scale, offset]).all():
            raise ValueError("Linked parameter paths must be unique and their scale/offset finite")
        used.add(path)
        normalized_links.append({"path": path, "scale": scale, "offset": offset})
    rows = []
    for index, driver_value in enumerate(values):
        run_jobs.progress(index, len(values), "Linked observable sweep")
        actual = {driver: float(driver_value)}
        try:
            current = set_parameter(model, driver, float(driver_value))
            for link in normalized_links:
                linked_value = link["scale"]*float(driver_value)+link["offset"]
                current = set_parameter(current, link["path"], linked_value)
                actual[link["path"]] = linked_value
            result = solve_stack(current)
            value = observable_value(result, spec)
            validity = "Not checked"
            convergence = None
            if validate_each:
                convergence = stack_convergence(current, tolerance=tolerance,
                    observables=[spec | {"tolerance": tolerance}])
                validity = "Validated" if convergence["converged"] and convergence["physical_balance_ok"] else "Not converged"
            submitted = asdict(current)
            rows.append({"driver_value": float(driver_value), "parameters": actual,
                         "value": value, "R": result["R"], "T": result["T"],
                         "A": result["A"], "validity": validity,
                         "convergence": convergence,
                         "submitted_model_sha256": settings_fingerprint(submitted)})
        except (ValueError, np.linalg.LinAlgError) as exc:
            rows.append({"driver_value": float(driver_value), "parameters": actual,
                         "value": None, "validity": "Failed", "error": str(exc)})
    solved = [row for row in rows if row["value"] is not None]
    eligible = [row for row in solved if not validate_each or row["validity"] == "Validated"]
    if not eligible:
        raise ValueError("No eligible linked-sweep points were computed; inspect failed or unconverged points")
    best = (max if goal == "max" else min)(eligible, key=lambda row: row["value"])
    return {"driver": driver, "driver_label": parameter_label(driver),
            "links": normalized_links, "observable": spec, "goal": goal,
            "rows": rows, "best": best, "validate_each": validate_each,
            "tolerance": tolerance, "submitted_model": asdict(model),
            "submitted_model_sha256": settings_fingerprint(asdict(model)),
            "status_counts": {status: sum(row["validity"] == status for row in rows)
                              for status in ("Validated", "Not converged", "Not checked", "Failed")}}
