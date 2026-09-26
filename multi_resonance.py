"""Joint far-field line-shape fitting for spectra containing several features.

This module fits sampled power spectra.  Its components are phenomenological
far-field features and are deliberately not called poles or quasinormal modes.
"""
from __future__ import annotations

from dataclasses import asdict
from itertools import product

import numpy as np
from scipy.optimize import curve_fit, linear_sum_assignment
from scipy.signal import find_peaks, peak_prominences

from stack import StackModel, solve_stack, stack_convergence, vertical_field
from sweep import set_parameter, sweep_values


def _line(x, center, width, amplitude, q, fano):
    e = 2*(x-center)/width
    if fano:
        # Subtract the asymptotic value so the shared polynomial is the
        # background and each resonant contribution tends to zero.
        return amplitude*((q+e)**2/(1+e**2)-1)
    return amplitude/(1+e**2)


def _candidate_indices(x, y, extremum, maximum):
    scale = max(float(np.ptp(y)), 1e-12)
    distance = max(2, len(x)//80)
    found = []
    signs = (1, -1) if extremum == "both" else ((1,) if extremum == "max" else (-1,))
    for sign in signs:
        indices, _ = find_peaks(sign*y, prominence=max(scale*.01, 1e-8), distance=distance)
        if indices.size:
            prominence = peak_prominences(sign*y, indices)[0]
            found.extend((int(i), float(p), "maximum" if sign == 1 else "minimum")
                         for i, p in zip(indices, prominence))
    found.sort(key=lambda row: row[1], reverse=True)
    selected = []
    for row in found:
        if all(abs(row[0]-old[0]) > 1 for old in selected):
            selected.append(row)
        if len(selected) >= maximum:
            break
    return sorted(selected, key=lambda row: x[row[0]])


def fit_multi_resonance_spectrum(wavelength_um, values, maximum_resonances=4,
                                 extremum="both", background_degree=1):
    """Select and jointly fit one or more Lorentzian/Fano far-field features."""
    x = np.asarray(wavelength_um, dtype=float)
    y = np.asarray(values, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 15 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Use at least 15 finite wavelength/value samples")
    if not np.all(np.diff(x) > 0):
        raise ValueError("Wavelength samples must be strictly increasing")
    if extremum not in ("max", "min", "both") or not 1 <= maximum_resonances <= 6:
        raise ValueError("Choose max, min, or both and 1–6 resonances")
    if background_degree not in (0, 1, 2):
        raise ValueError("Background degree must be 0, 1, or 2")
    candidates = _candidate_indices(x, y, extremum, maximum_resonances)
    if not candidates:
        raise ValueError("No feature passed the global prominence screen")
    span, dx = float(x[-1]-x[0]), float(np.median(np.diff(x)))
    xmid = float(np.mean(x))
    models = []
    for count in range(1, len(candidates)+1):
        chosen = sorted(candidates, key=lambda row: row[1], reverse=True)[:count]
        chosen.sort(key=lambda row: x[row[0]])
        for shapes in product((False, True), repeat=count):
            nbg = background_degree+1
            p0 = np.polyfit(x-xmid, y, background_degree)[::-1].tolist()
            lower = [-2, *([-1e5]*background_degree)]
            upper = [2, *([1e5]*background_degree)]
            for (index, _, kind), fano in zip(chosen, shapes):
                local = y[index]-float(np.median(y))
                amplitude = local if abs(local) > 1e-6 else (.05 if kind == "maximum" else -.05)
                p0.extend([x[index], max(3*dx, span/50), amplitude, 1.0])
                lower.extend([x[0], max(dx*.25, span/1e6), -5, -50 if fano else -.001])
                upper.extend([x[-1], span, 5, 50 if fano else .001])

            def evaluate(xx, *parameters):
                result = sum(parameters[j]*(xx-xmid)**j for j in range(nbg))
                cursor = nbg
                for fano in shapes:
                    center, width, amplitude, q = parameters[cursor:cursor+4]
                    result = result+_line(xx, center, width, amplitude, q, fano)
                    cursor += 4
                return result
            try:
                popt, pcov = curve_fit(evaluate, x, y, p0=p0, bounds=(lower, upper),
                                       maxfev=100000, x_scale="jac")
            except (RuntimeError, ValueError, FloatingPointError):
                continue
            predicted = evaluate(x, *popt)
            residual = y-predicted
            rss = max(float(np.dot(residual, residual)), 1e-30)
            p = len(popt)
            bic = len(x)*np.log(rss/len(x))+p*np.log(len(x))
            sigma = np.sqrt(np.maximum(0, np.diag(pcov)))
            correlation = pcov/np.maximum(np.outer(sigma, sigma), 1e-300)
            components = []
            for j, fano in enumerate(shapes):
                base = nbg+4*j
                center, width, amplitude, q = popt[base:base+4]
                center_sigma, width_sigma = sigma[base], sigma[base+1]
                components.append({"index": j+1, "line_shape": "Fano" if fano else "Lorentzian",
                    "center_um": float(center), "linewidth_um": float(abs(width)),
                    "loaded_Q": float(abs(center/width)), "amplitude": float(amplitude),
                    "fano_q": float(q) if fano else None,
                    "center_sigma_um": float(center_sigma),
                    "linewidth_sigma_um": float(width_sigma), "resolved": True,
                    "warnings": []})
            warnings = []
            for i, component in enumerate(components):
                if component["linewidth_sigma_um"] > .5*component["linewidth_um"]:
                    component["resolved"] = False
                    component["warnings"].append("linewidth uncertainty exceeds 50%")
                for j in range(i):
                    other = components[j]
                    separation = abs(component["center_um"]-other["center_um"])
                    if separation < .5*(component["linewidth_um"]+other["linewidth_um"]):
                        component["resolved"] = other["resolved"] = False
                        message = f"overlaps component {j+1} within the mean linewidth"
                        component["warnings"].append(message)
                        other["warnings"].append(f"overlaps component {i+1} within the mean linewidth")
                    ci, cj = nbg+4*i, nbg+4*j
                    block = np.max(np.abs(correlation[ci:ci+2, cj:cj+2]))
                    if block > .95:
                        component["resolved"] = other["resolved"] = False
                        component["warnings"].append(f"strong parameter correlation with component {j+1}")
                        other["warnings"].append(f"strong parameter correlation with component {i+1}")
            if any(not c["resolved"] for c in components):
                warnings.append("One or more fitted components are not individually identifiable")
            models.append({"resonance_count": count, "shapes": ["Fano" if s else "Lorentzian" for s in shapes],
                "bic": float(bic), "rmse": float(np.sqrt(rss/len(x))), "parameters": popt.tolist(),
                "parameter_sigma": sigma.tolist(), "correlation": correlation.tolist(),
                "components": components, "predicted": predicted.tolist(),
                "residual": residual.tolist(), "warnings": warnings})
    if not models:
        raise ValueError("No joint resonance model converged")
    models.sort(key=lambda item: item["bic"])
    best = models[0]
    delta = float(models[1]["bic"]-best["bic"]) if len(models) > 1 else None
    if delta is not None and delta < 2:
        best["warnings"].append("Competing line-shape/count models have ΔBIC < 2")
    return {"wavelength_um": x.tolist(), "values": y.tolist(),
            "candidates": [{"wavelength_um": float(x[i]), "prominence": p, "kind": kind}
                           for i, p, kind in candidates],
            "best_model": best,
            "model_comparison": [{"resonance_count": m["resonance_count"], "shapes": m["shapes"],
                                  "bic": m["bic"], "delta_bic": float(m["bic"]-best["bic"]),
                                  "rmse": m["rmse"]} for m in models],
            "evidence_type": "joint phenomenological far-field line-shape fit",
            "not_provided": ["complex-frequency poles", "quasinormal modes", "radiative/absorptive Q decomposition"]}


def adaptive_multi_resonance(model: StackModel, start, stop, quantity="R", extremum="both",
                             global_points=101, refinement_points=31, rounds=2,
                             maximum_resonances=4, background_degree=1):
    """Global scan, candidate-wise refinement, then a joint line-shape fit."""
    if quantity not in ("R", "T", "A") or not 31 <= global_points <= 401:
        raise ValueError("Choose R, T, or A and 31–401 global points")
    if not 15 <= refinement_points <= 101 or not 0 <= rounds <= 4 or not 0 < start < stop:
        raise ValueError("Use increasing wavelengths, 15–101 refinement points, and 0–4 rounds")
    cache = {}
    def sample(wavelength):
        key = round(float(wavelength), 14)
        if key not in cache:
            current = StackModel(**{**asdict(model), "layers": model.layers, "wavelength_um": float(wavelength)})
            result = solve_stack(current)
            cache[key] = {"wavelength_um": float(wavelength), **{k: float(result[k]) for k in ("R", "T", "A")}}
        return cache[key]
    for wavelength in np.linspace(start, stop, global_points):
        sample(wavelength)
    for _ in range(rounds):
        rows = sorted(cache.values(), key=lambda row: row["wavelength_um"])
        x = np.asarray([r["wavelength_um"] for r in rows]); y = np.asarray([r[quantity] for r in rows])
        candidates = _candidate_indices(x, y, extremum, maximum_resonances)
        if not candidates:
            break
        for index, _, _kind in candidates:
            half = max((stop-start)/(global_points-1)*3/(2**_), np.median(np.diff(x))*2)
            lo, hi = max(start, x[index]-half), min(stop, x[index]+half)
            for wavelength in np.linspace(lo, hi, refinement_points):
                sample(wavelength)
    rows = sorted(cache.values(), key=lambda row: row["wavelength_um"])
    fitted = fit_multi_resonance_spectrum([r["wavelength_um"] for r in rows],
        [r[quantity] for r in rows], maximum_resonances, extremum, background_degree)
    fitted["rows"] = rows
    fitted["quantity"] = quantity
    fitted["convergence"] = []
    for component in fitted["best_model"]["components"]:
        current = StackModel(**{**asdict(model), "layers": model.layers,
                               "wavelength_um": component["center_um"]})
        check = stack_convergence(current)
        field = vertical_field(current, "xz", .5, 21, 7)
        signature = np.sqrt(np.maximum(np.asarray(field["E2"], dtype=float), 0)).ravel()
        signature /= max(float(np.linalg.norm(signature)), 1e-30)
        component["field_signature"] = signature.tolist()
        component["field_signature_scope"] = "normalized coarse x-z |E| amplitude; branch-tracking diagnostic"
        fitted["convergence"].append({"center_um": component["center_um"], "result": check})
    return fitted


def track_resonance_branches(steps, field_overlaps=None):
    """Match fitted features across a sweep using wavelength and optional field overlap."""
    if not steps:
        return {"branches": [], "warnings": []}
    branches = [[dict(component, step=0)] for component in steps[0]]
    warnings = []
    for step_index, current in enumerate(steps[1:], 1):
        previous = [branch[-1] for branch in branches]
        if not previous or not current:
            warnings.append(f"step {step_index}: no one-to-one branch match")
            continue
        cost = np.zeros((len(previous), len(current)))
        for i, old in enumerate(previous):
            for j, new in enumerate(current):
                wavelength_scale = max(old["linewidth_um"], new["linewidth_um"], 1e-12)
                cost[i, j] = abs(old["center_um"]-new["center_um"])/wavelength_scale
                if field_overlaps is not None:
                    overlap = float(field_overlaps[step_index-1][i][j])
                    cost[i, j] += 2*(1-np.clip(overlap, 0, 1))
                elif old.get("field_signature") is not None and new.get("field_signature") is not None:
                    a, b = np.asarray(old["field_signature"]), np.asarray(new["field_signature"])
                    overlap = float(abs(np.vdot(a, b))**2/(max(np.vdot(a, a).real*np.vdot(b, b).real, 1e-30)))
                    cost[i, j] += 2*(1-np.clip(overlap, 0, 1))
        row, col = linear_sum_assignment(cost)
        for i, j in zip(row, col):
            item = dict(current[j], step=step_index, tracking_cost=float(cost[i, j]))
            if cost[i, j] > 2:
                item.setdefault("warnings", []).append("weak branch match")
                warnings.append(f"step {step_index}: branch {i+1} has weak continuity")
            branches[i].append(item)
        for j in set(range(len(current)))-set(col):
            branches.append([dict(current[j], step=step_index, warnings=["new/unmatched branch"])])
    return {"branches": branches, "warnings": warnings,
            "method": "Hungarian assignment using linewidth-normalized wavelength distance and optional normalized field overlap"}


def multi_resonance_sweep(model: StackModel, parameter: str, start_value: float,
                          stop_value: float, parameter_points: int, wavelength_start: float,
                          wavelength_stop: float, quantity="R", extremum="both",
                          global_points=61, refinement_points=21, rounds=1,
                          maximum_resonances=4, background_degree=1):
    """Fit all visible features at each sweep step and assign continuous branches."""
    if not 2 <= parameter_points <= 15:
        raise ValueError("Use 2–15 branch-tracking parameter points")
    values = sweep_values(start_value, stop_value, parameter_points)
    steps, results = [], []
    for value in values:
        current = set_parameter(model, parameter, float(value))
        result = adaptive_multi_resonance(current, wavelength_start, wavelength_stop,
            quantity, extremum, global_points, refinement_points, rounds,
            maximum_resonances, background_degree)
        components = result["best_model"]["components"]
        steps.append(components)
        results.append({"parameter": float(value), "components": components,
                        "model_comparison": result["model_comparison"],
                        "warnings": result["best_model"]["warnings"]})
    tracking = track_resonance_branches(steps)
    return {"parameter": parameter, "parameter_values": values.tolist(),
            "steps": results, **tracking,
            "evidence_type": "far-field fitted branches with coarse intensity-overlap continuity",
            "warning": "Branch assignment is diagnostic. Avoid claiming an avoided crossing, BIC, or modal identity without converged eigenmode or pole fields."}
