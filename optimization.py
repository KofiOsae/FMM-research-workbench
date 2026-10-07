"""Bounded multi-objective geometry optimization for the finite-stack solver."""

from __future__ import annotations

from dataclasses import asdict
import numpy as np
from scipy.optimize import differential_evolution, LinearConstraint
from run_jobs import progress

from stack import StackModel, solve_stack
from sweep import set_parameter
from observables import (design_metric_spec, design_metric_value,
                         evaluate_named_observables,
                         normalize_observable_definitions)
from research_validation import validation_report


def optimize_geometry(model: StackModel, variables: list[dict], objectives: list[dict],
                      linear_constraints: list[dict] | None = None,
                      generations: int = 4, population: int = 5,
                      seed: int = 12345, polish: bool = True,
                      robust_samples: int = 1, robust_weight: float = .25,
                      convergence_aware: bool = False,
                      convergence_tolerance: float = .01,
                      observable_definitions: list[dict] | None = None,
                      optical_constraints: list[dict] | None = None,
                      constraint_tolerance: float = 1e-6) -> dict:
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
    if not 0 <= constraint_tolerance <= .05:
        raise ValueError("Constraint feasibility tolerance must be from 0 to 0.05")
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
    definitions = normalize_observable_definitions(observable_definitions)
    normalized_objectives = []
    for item in objectives:
        wavelength = float(item["wavelength_um"]); spec = design_metric_spec(item, definitions)
        goal = str(item.get("goal", "max")); weight = float(item.get("weight", 1))
        target = float(item.get("target", 0))
        if wavelength <= 0 or goal not in ("max", "min", "target") or weight <= 0 or not np.isfinite([wavelength, weight, target]).all():
            raise ValueError("Each objective needs positive wavelength and weight, supported quantity, and max/min/target goal")
        normalized_objectives.append({"wavelength_um": wavelength, **spec,
                                      "goal": goal, "weight": weight, "target": target})
    normalized_optical_constraints = []
    if len(optical_constraints or []) > 12:
        raise ValueError("Use no more than 12 optical constraints")
    for item in optical_constraints or []:
        wavelength = float(item["wavelength_um"])
        lower = float(item.get("lower", -np.inf)) if item.get("lower") not in (None, "") else -np.inf
        upper = float(item.get("upper", np.inf)) if item.get("upper") not in (None, "") else np.inf
        spec = design_metric_spec(item, definitions)
        if (wavelength <= 0 or lower > upper or
                not np.isfinite(wavelength) or (not np.isfinite(lower) and not np.isfinite(upper))):
            raise ValueError("Each optical constraint needs a positive wavelength and at least one valid bound")
        normalized_optical_constraints.append({"wavelength_um": wavelength, **spec,
            "lower": None if not np.isfinite(lower) else lower,
            "upper": None if not np.isfinite(upper) else upper,
            "label": str(item.get("label") or spec["label"])[:120]})
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

    def metric_rows(current):
        """Solve each requested wavelength once and evaluate all named metrics."""
        wavelengths = sorted({item["wavelength_um"] for item in
                              normalized_objectives+normalized_optical_constraints})
        solved = {}
        for wavelength in wavelengths:
            point = StackModel(**{**asdict(current), "layers": current.layers,
                "wavelength_um": wavelength})
            result_point = solve_stack(point)
            solved[wavelength] = (point, result_point,
                                  evaluate_named_observables(result_point, definitions))
        terms, total = [], 0.0
        for objective in normalized_objectives:
            point, result_point, named = solved[objective["wavelength_um"]]
            value = design_metric_value(result_point, objective, definitions, named)
            raw = (-value if objective["goal"] == "max" else value if
                   objective["goal"] == "min" else (value-objective["target"])**2)
            contribution = objective["weight"]*raw
            convergence_change = None
            if convergence_aware:
                high_budget = min(201, point.order_budget+24)
                if high_budget == point.order_budget:
                    raise ValueError("Convergence-aware optimization needs an order budget below 201")
                high = StackModel(**{**asdict(point), "layers": point.layers,
                                     "order_budget": high_budget})
                high_result = solve_stack(high)
                high_named = evaluate_named_observables(high_result, definitions)
                high_value = design_metric_value(high_result, objective, definitions, high_named)
                convergence_change = abs(high_value-value)
                contribution += 1e3*max(0.0, convergence_change-convergence_tolerance)
            total += contribution
            terms.append({**objective, "value": value, "loss_contribution": contribution,
                          "convergence_change": convergence_change})
        constraint_terms, constraint_penalty = [], 0.0
        for constraint in normalized_optical_constraints:
            _, result_point, named = solved[constraint["wavelength_um"]]
            value = design_metric_value(result_point, constraint, definitions, named)
            below = (max(0.0, constraint["lower"]-value)
                     if constraint["lower"] is not None else 0.0)
            above = (max(0.0, value-constraint["upper"])
                     if constraint["upper"] is not None else 0.0)
            violation = max(below, above)
            effective = max(0.0, violation-constraint_tolerance)
            penalty = 1e4*effective+1e6*effective*effective
            constraint_penalty += penalty
            constraint_terms.append({**constraint, "value": value,
                "violation": violation, "feasible": violation <= constraint_tolerance,
                "penalty": penalty})
        return terms, constraint_terms, float(total+constraint_penalty), float(constraint_penalty)

    def evaluate(x, record=True):
        nonlocal failures
        sample_losses, sample_terms, sample_constraints, sample_penalties = [], [], [], []
        for offset in robust_offsets:
            trial = np.clip(np.asarray(x)+offset, np.asarray(bounds)[:,0], np.asarray(bounds)[:,1])
            current, terms, constraint_terms, total, penalty = build(trial), [], [], 0.0, 0.0
            try:
                terms, constraint_terms, total, penalty = metric_rows(current)
            except (ValueError, np.linalg.LinAlgError):
                failures += 1; total = 1e9; terms = []; constraint_terms = []
            sample_losses.append(float(total)); sample_terms.append(terms)
            sample_constraints.append(constraint_terms); sample_penalties.append(penalty)
        total = float(np.mean(sample_losses)+robust_weight*np.std(sample_losses))
        terms = sample_terms[0]
        nominal_constraints = sample_constraints[0]
        worst_violation = max((row["violation"] for rows in sample_constraints for row in rows), default=0.0)
        if record: evaluations.append({"parameters": {p: float(v) for p,v in zip(paths,x)},
                                       "loss": total, "objectives": terms,
                                       "constraints": nominal_constraints,
                                       "feasible": bool(nominal_constraints and all(r["feasible"] for r in nominal_constraints)) if normalized_optical_constraints else True,
                                       "robust_feasible": all(r["feasible"] for rows in sample_constraints for r in rows),
                                       "maximum_constraint_violation": worst_violation,
                                       "constraint_penalty": float(np.mean(sample_penalties)),
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
    dashboard = []
    for wavelength in sorted({item["wavelength_um"] for item in
                              normalized_objectives+normalized_optical_constraints}):
        point = StackModel(**{**asdict(best_model), "layers": best_model.layers,
            "wavelength_um": wavelength})
        scattering = solve_stack(point)
        named = evaluate_named_observables(scattering, definitions)
        dashboard.append({"wavelength_um": wavelength,
            "R": scattering["R"], "T": scattering["T"], "A": scattering["A"],
            "R0": scattering["R0"], "T0": scattering["T0"],
            "observables": named,
            "propagating_orders": [{"m": row["m"], "n": row["n"],
                "R": row["R"], "T": row["T"], "angle_reflected_deg": row.get("angle_reflected_deg"),
                "angle_transmitted_deg": row.get("angle_transmitted_deg")}
                for row in scattering.get("orders", [])
                if row["R"] > 1e-12 or row["T"] > 1e-12]})
    checks = []
    requested_items = normalized_objectives+normalized_optical_constraints
    for wavelength in sorted({item["wavelength_um"] for item in requested_items}):
        point = StackModel(**{**asdict(best_model), "layers": best_model.layers,
            "wavelength_um": wavelength})
        requested = [{"quantity": "expression", "expression": item["expression"],
                      "label": item["label"], "tolerance": convergence_tolerance}
                     for item in requested_items if item["wavelength_um"] == wavelength]
        trust = validation_report(point, requested, tolerance=convergence_tolerance,
                                  adaptive=True)
        checks.append({"wavelength_um": wavelength, "status": trust["status"],
            "checks": trust["checks"], "fourier_change": trust["fourier"]["max_change"],
            "grid_change": (trust["geometry_grid"] or {}).get("max_change"),
            "submitted_model_sha256": trust["submitted_model_sha256"],
            "requested_observables": trust["fourier"].get("requested_observables", []),
            "order_map": trust["order_map"]})
    feasible = bool(best_evaluation["feasible"])
    trust_validated = all(row["status"] == "Validated" for row in checks)
    return {"success": bool(result.success), "message": str(result.message),
        "best_loss": final, "best_parameters": {p: float(v) for p,v in zip(paths,result.x)},
        "best_model": asdict(best_model), "observables": definitions,
        "objectives": normalized_objectives, "optical_constraints": normalized_optical_constraints,
        "linear_constraints": constraint_records, "objective_results": best_evaluation["objectives"],
        "constraint_results": best_evaluation["constraints"], "feasible": feasible,
        "design_dashboard": dashboard,
        "robust_feasible": best_evaluation["robust_feasible"],
        "maximum_constraint_violation": best_evaluation["maximum_constraint_violation"],
        "evaluations": evaluations, "evaluation_count": int(result.nfev), "failed_evaluations": failures,
        "history": history, "seed": seed, "generations": generations, "population": population,
        "polished": bool(polish), "validation": checks, "trust_validation": checks,
        "robustness":{"samples":robust_samples,"weight":robust_weight,
                      "uncertainty_sigma":{path:value for path,value in zip(paths,uncertainty_sigma)},
                      "mean_loss":best_evaluation["robust_mean_loss"],
                      "std_loss":best_evaluation["robust_std_loss"],
                      "worst_loss":best_evaluation["robust_worst_loss"]},
        "convergence_aware_candidates": bool(convergence_aware),
        "convergence_tolerance": convergence_tolerance,
        "validated": bool(feasible and trust_validated),
        "constraint_tolerance": constraint_tolerance,
        "warning": "This bounded stochastic search does not prove a global optimum. Feasibility and Trust validation are reported separately. Robust samples use fixed Gaussian fabrication perturbations clipped to the design bounds. Repeat with multiple seeds, justified distributions, and reserved validation wavelengths."}
