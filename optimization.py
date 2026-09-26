"""Bounded multi-objective geometry optimization for the finite-stack solver."""

from __future__ import annotations

from dataclasses import asdict
import numpy as np
from scipy.optimize import differential_evolution, LinearConstraint

from stack import StackModel, solve_stack, stack_convergence, grid_convergence
from sweep import set_parameter


def optimize_geometry(model: StackModel, variables: list[dict], objectives: list[dict],
                      linear_constraints: list[dict] | None = None,
                      generations: int = 4, population: int = 5,
                      seed: int = 12345, polish: bool = True) -> dict:
    model.validate()
    if not 1 <= len(variables) <= 5:
        raise ValueError("Use one to five design variables")
    if not 1 <= len(objectives) <= 8:
        raise ValueError("Use one to eight optical objectives")
    if not 1 <= generations <= 20 or not 4 <= population <= 12:
        raise ValueError("Use 1–20 generations and population multiplier 4–12")
    paths, bounds = [], []
    for item in variables:
        path, lower, upper = str(item["path"]), float(item["lower"]), float(item["upper"])
        if path in paths or not np.isfinite([lower, upper]).all() or lower >= upper:
            raise ValueError("Design-variable paths must be unique with increasing finite bounds")
        set_parameter(model, path, (lower+upper)/2)
        paths.append(path); bounds.append((lower, upper))
    normalized_objectives = []
    for item in objectives:
        wavelength = float(item["wavelength_um"]); quantity = str(item["quantity"])
        goal = str(item.get("goal", "max")); weight = float(item.get("weight", 1))
        target = float(item.get("target", 0))
        if wavelength <= 0 or quantity not in ("R", "T", "A", "R0", "T0") or goal not in ("max", "min", "target") or weight <= 0 or not np.isfinite([wavelength, weight, target]).all():
            raise ValueError("Each objective needs positive wavelength and weight, supported quantity, and max/min/target goal")
        normalized_objectives.append({"wavelength_um": wavelength, "quantity": quantity,
                                      "goal": goal, "weight": weight, "target": target})
    scipy_constraints = []
    constraint_records = []
    for item in linear_constraints or []:
        coefficients = np.asarray(item["coefficients"], dtype=float)
        lower, upper = float(item.get("lower", -np.inf)), float(item.get("upper", np.inf))
        if coefficients.shape != (len(paths),) or np.isnan(coefficients).any() or lower > upper:
            raise ValueError("Every linear constraint needs one coefficient per variable and lower ≤ upper")
        scipy_constraints.append(LinearConstraint(coefficients, lower, upper))
        constraint_records.append({"coefficients": coefficients.tolist(),
            "lower": None if not np.isfinite(lower) else lower,
            "upper": None if not np.isfinite(upper) else upper})
    evaluations, failures = [], 0
    def build(x):
        current = model
        for path, value in zip(paths, x): current = set_parameter(current, path, float(value))
        return current
    def evaluate(x, record=True):
        nonlocal failures
        current, terms, total = build(x), [], 0.0
        try:
            for objective in normalized_objectives:
                point = StackModel(**{**asdict(current), "layers": current.layers,
                    "wavelength_um": objective["wavelength_um"]})
                value = float(solve_stack(point)[objective["quantity"]])
                raw = -value if objective["goal"] == "max" else value if objective["goal"] == "min" else (value-objective["target"])**2
                contribution = objective["weight"]*raw
                total += contribution
                terms.append({**objective, "value": value, "loss_contribution": contribution})
        except (ValueError, np.linalg.LinAlgError):
            failures += 1; total = 1e6; terms = []
        if record: evaluations.append({"parameters": {p: float(v) for p,v in zip(paths,x)},
                                       "loss": float(total), "objectives": terms})
        return float(total)
    history = []
    def callback(intermediate_result):
        history.append({"generation": len(history)+1, "best_loss": float(intermediate_result.fun),
                        "parameters": {p: float(v) for p,v in zip(paths, intermediate_result.x)}})
    result = differential_evolution(evaluate, bounds, constraints=tuple(scipy_constraints),
        maxiter=generations, popsize=population, seed=seed, polish=polish,
        callback=callback, updating="immediate", workers=1, tol=1e-4)
    best_model = build(result.x)
    final = evaluate(result.x, record=True)
    best_evaluation = evaluations[-1]
    checks = []
    for objective in normalized_objectives:
        point = StackModel(**{**asdict(best_model), "layers": best_model.layers,
            "wavelength_um": objective["wavelength_um"]})
        order = stack_convergence(point); grid = grid_convergence(point)
        checks.append({"wavelength_um": objective["wavelength_um"],
            "fourier_converged": bool(order["converged"] and order["physical_balance_ok"]),
            "grid_converged": bool(grid["converged"]), "fourier_change": order["max_change"],
            "grid_change": grid["max_change"]})
    return {"success": bool(result.success), "message": str(result.message),
        "best_loss": final, "best_parameters": {p: float(v) for p,v in zip(paths,result.x)},
        "best_model": asdict(best_model), "objectives": normalized_objectives,
        "linear_constraints": constraint_records, "objective_results": best_evaluation["objectives"],
        "evaluations": evaluations, "evaluation_count": int(result.nfev), "failed_evaluations": failures,
        "history": history, "seed": seed, "generations": generations, "population": population,
        "polished": bool(polish), "validation": checks,
        "validated": all(row["fourier_converged"] and row["grid_converged"] for row in checks),
        "warning": "This bounded stochastic search does not prove a global optimum. Repeat with multiple seeds, wider justified bounds, tolerance analysis, and reserved validation wavelengths."}
