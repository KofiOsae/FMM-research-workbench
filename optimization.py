"""Bounded multi-objective geometry optimization for the finite-stack solver."""

from __future__ import annotations

from dataclasses import asdict
import numpy as np
from scipy.optimize import differential_evolution, LinearConstraint
from run_jobs import progress

from stack import StackModel, solve_stack, stack_convergence, grid_convergence
from sweep import set_parameter
from observables import observable_spec, observable_value


def optimize_geometry(model: StackModel, variables: list[dict], objectives: list[dict],
                      linear_constraints: list[dict] | None = None,
                      generations: int = 4, population: int = 5,
                      seed: int = 12345, polish: bool = True,
                      robust_samples: int = 1, robust_weight: float = .25,
                      convergence_aware: bool = False,
                      convergence_tolerance: float = .01) -> dict:
    model.validate()
    if not 1 <= len(variables) <= 5:
        raise ValueError("Use one to five design variables")
    if not 1 <= len(objectives) <= 8:
        raise ValueError("Use one to eight optical objectives")
    if not 1 <= generations <= 20 or not 4 <= population <= 12:
        raise ValueError("Use 1–20 generations and population multiplier 4–12")
    if not 1 <= robust_samples <= 12 or not 0 <= robust_weight <= 2:
        raise ValueError("Use 1–12 robust samples and robustness weight from 0 to 2")
    if not 1e-8 <= convergence_tolerance <= .25:
        raise ValueError("Convergence tolerance must be from 1e-8 to 0.25")
    paths, bounds, uncertainty_sigma = [], [], []
    for item in variables:
        path, lower, upper = str(item["path"]), float(item["lower"]), float(item["upper"])
        if path in paths or not np.isfinite([lower, upper]).all() or lower >= upper:
            raise ValueError("Design-variable paths must be unique with increasing finite bounds")
        set_parameter(model, path, (lower+upper)/2)
        sigma = float(item.get("uncertainty_sigma", 0))
        if not np.isfinite(sigma) or sigma < 0:
            raise ValueError("Fabrication uncertainty sigma must be finite and nonnegative")
        paths.append(path); bounds.append((lower, upper)); uncertainty_sigma.append(sigma)
    normalized_objectives = []
    for item in objectives:
        wavelength = float(item["wavelength_um"]); spec = observable_spec(item)
        goal = str(item.get("goal", "max")); weight = float(item.get("weight", 1))
        target = float(item.get("target", 0))
        if wavelength <= 0 or goal not in ("max", "min", "target") or weight <= 0 or not np.isfinite([wavelength, weight, target]).all():
            raise ValueError("Each objective needs positive wavelength and weight, supported quantity, and max/min/target goal")
        normalized_objectives.append({"wavelength_um": wavelength, **spec,
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
    rng = np.random.default_rng(seed+7919)
    robust_offsets = np.zeros((robust_samples, len(paths)))
    if robust_samples > 1:
        robust_offsets[1:] = rng.normal(size=(robust_samples-1, len(paths)))*np.asarray(uncertainty_sigma)
    def build(x):
        current = model
        for path, value in zip(paths, x): current = set_parameter(current, path, float(value))
        return current
    def evaluate(x, record=True):
        nonlocal failures
        sample_losses, sample_terms = [], []
        for offset in robust_offsets:
            trial = np.clip(np.asarray(x)+offset, np.asarray(bounds)[:,0], np.asarray(bounds)[:,1])
            current, terms, total = build(trial), [], 0.0
            try:
                for objective in normalized_objectives:
                    point = StackModel(**{**asdict(current), "layers": current.layers,
                        "wavelength_um": objective["wavelength_um"]})
                    result_point = solve_stack(point)
                    value = observable_value(result_point, objective)
                    raw = -value if objective["goal"] == "max" else value if objective["goal"] == "min" else (value-objective["target"])**2
                    contribution = objective["weight"]*raw
                    convergence_change = None
                    if convergence_aware:
                        high_budget = min(201, point.order_budget+24)
                        if high_budget == point.order_budget:
                            raise ValueError("Convergence-aware optimization needs an order budget below 201")
                        high = StackModel(**{**asdict(point), "layers": point.layers,
                                             "order_budget": high_budget})
                        high_value = observable_value(solve_stack(high), objective)
                        convergence_change = abs(high_value-value)
                        contribution += 1e3*max(0.0, convergence_change-convergence_tolerance)
                    total += contribution
                    terms.append({**objective, "value": value, "loss_contribution": contribution,
                                  "convergence_change": convergence_change})
            except (ValueError, np.linalg.LinAlgError):
                failures += 1; total = 1e6; terms = []
            sample_losses.append(float(total)); sample_terms.append(terms)
        total = float(np.mean(sample_losses)+robust_weight*np.std(sample_losses))
        terms = sample_terms[0]
        if record: evaluations.append({"parameters": {p: float(v) for p,v in zip(paths,x)},
                                       "loss": total, "objectives": terms,
                                       "robust_mean_loss": float(np.mean(sample_losses)),
                                       "robust_std_loss": float(np.std(sample_losses)),
                                       "robust_worst_loss": float(np.max(sample_losses))})
        progress(len(evaluations), max(1, generations*population*len(paths)), "Geometry optimization")
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
        requested = [objective | {"tolerance": convergence_tolerance}]
        order = stack_convergence(point, tolerance=convergence_tolerance,
                                  observables=requested)
        grid = grid_convergence(point, tolerance=convergence_tolerance,
                                observables=requested)
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
        "robustness":{"samples":robust_samples,"weight":robust_weight,
                      "uncertainty_sigma":{path:value for path,value in zip(paths,uncertainty_sigma)},
                      "mean_loss":best_evaluation["robust_mean_loss"],
                      "std_loss":best_evaluation["robust_std_loss"],
                      "worst_loss":best_evaluation["robust_worst_loss"]},
        "convergence_aware_candidates": bool(convergence_aware),
        "convergence_tolerance": convergence_tolerance,
        "validated": all(row["fourier_converged"] and row["grid_converged"] for row in checks),
        "warning": "This bounded stochastic search does not prove a global optimum. Robust samples use fixed Gaussian fabrication perturbations clipped to the design bounds. Repeat with multiple seeds, justified distributions, and reserved validation wavelengths."}
