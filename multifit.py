"""Bounded multi-parameter fitting with covariance and validation diagnostics."""

from __future__ import annotations

from dataclasses import asdict, replace

from run_jobs import progress

import numpy as np
from scipy.optimize import least_squares, differential_evolution

from experiment import parse_measurement
from stack import StackModel, solve_stack
from sweep import set_parameter, parameter_label
from tmm import solve_tmm


ALLOWED = {"theta_deg", "incident_n", "exit_n", "period_x_um", "period_y_um"}
NUISANCE_KEYS = {"wavelength_offset_um", "scale", "background", "fwhm_um"}
NUISANCE_LABELS = {"wavelength_offset_um": "wavelength offset (µm)",
                   "scale": "intensity scale", "background": "intensity background",
                   "fwhm_um": "instrument FWHM (µm)"}


def _label(path: str) -> str:
    parts = path.split(".")
    if len(parts) == 3 and parts[0] == "dataset":
        return f"Dataset {int(parts[1])+1} {NUISANCE_LABELS[parts[2]]}"
    return parameter_label(path)


def _allowed_path(path: str) -> bool:
    nuisance = path.split(".")
    return (len(nuisance) == 3 and nuisance[0] == "dataset" and
            nuisance[1].isdigit() and nuisance[2] in NUISANCE_KEYS) or path in ALLOWED or (path.startswith("layer.") and
        path.split(".")[-1] in {"thickness_um", "background_n", "feature_n",
                                "fill_x", "fill_y", "inner_radius", "offset_x", "offset_y"})


def _instrument_response(values: np.ndarray, wavelengths: np.ndarray, fwhm_um: float) -> np.ndarray:
    """Apply a unit-area Gaussian line-spread function on an arbitrary wavelength grid."""
    if fwhm_um <= 0:
        return values
    sigma = fwhm_um / (2*np.sqrt(2*np.log(2)))
    distance = wavelengths[:, None]-wavelengths[None, :]
    weights = np.exp(-.5*(distance/sigma)**2)
    weights /= np.maximum(weights.sum(axis=1, keepdims=True), 1e-300)
    return weights @ values


def _predict(model: StackModel, wavelengths, channels) -> np.ndarray:
    uniform = all(layer.kind == "uniform" for layer in model.layers)
    rows = []
    for wavelength in wavelengths:
        current = replace(model, wavelength_um=float(wavelength))
        result = solve_tmm(current) if uniform else solve_stack(current)
        rows.append([result[key] for key in channels])
    return np.asarray(rows)


def multi_parameter_fit(model: StackModel, datasets: list[dict], parameters: list[dict],
                        method: str = "local", validation_fraction: float = .2,
                        bootstrap: int = 0, seed: int = 12345) -> dict:
    model.validate()
    if not 1 <= len(datasets) <= 6:
        raise ValueError("Use one to six measured datasets")
    if not 1 <= len(parameters) <= 5:
        raise ValueError("Fit one to five parameters")
    if method not in ("local", "global_then_local"):
        raise ValueError("Fit method must be local or global_then_local")
    if not 0 <= validation_fraction <= .45 or not 0 <= bootstrap <= 200:
        raise ValueError("Validation fraction must be 0–0.45 and bootstrap 0–200")
    paths, initial, lower, upper = [], [], [], []
    for item in parameters:
        path = str(item["path"])
        if not _allowed_path(path):
            raise ValueError(f"Unsupported fit parameter: {path}")
        _label(path)
        lo, hi = float(item["lower"]), float(item["upper"])
        if not np.isfinite([lo, hi]).all() or lo >= hi:
            raise ValueError("Every fit parameter needs finite increasing bounds")
        if not path.startswith("dataset."):
            probe = set_parameter(model, path, (lo+hi)/2)
            del probe
        if path.startswith("dataset."):
            dataset_index, key = int(path.split(".")[1]), path.split(".")[2]
            if dataset_index >= len(datasets):
                raise ValueError(f"Nuisance parameter refers to missing dataset {dataset_index}")
            defaults = {"wavelength_offset_um": 0., "scale": 1., "background": 0., "fwhm_um": 0.}
            value = float(datasets[dataset_index].get(key, defaults[key]))
            if key == "fwhm_um" and lo < 0:
                raise ValueError("Instrument FWHM cannot be negative")
        elif path.startswith("layer."):
            _, raw, key = path.split(".")
            value = getattr(model.layers[int(raw)], key)
        else:
            value = getattr(model, path)
        if not lo <= value <= hi:
            raise ValueError(f"Initial {path} lies outside its bounds")
        paths.append(path); initial.append(float(value)); lower.append(lo); upper.append(hi)
    if len(set(paths)) != len(paths):
        raise ValueError("Fit parameters must be different")

    parsed = []
    for item in datasets:
        wavelengths, channels, measured = parse_measurement(item["csv"])
        polarization = str(item.get("polarization", model.polarization))
        if polarization not in ("s", "p", "unpolarized"):
            raise ValueError("Dataset polarization must be s, p, or unpolarized")
        order = np.arange(len(wavelengths))
        # Evenly distribute validation points across the spectrum.
        validation_count = int(round(validation_fraction*len(order)))
        validation_index = set(np.linspace(0, len(order)-1, validation_count,
                                           dtype=int).tolist()) if validation_count else set()
        train = np.asarray([i not in validation_index for i in order])
        parsed.append((wavelengths, channels, measured, polarization, train))

    def build(values):
        current = model
        for path, value in zip(paths, values):
            if not path.startswith("dataset."):
                current = set_parameter(current, path, float(value))
        return current

    def nuisance(values, dataset_index):
        result = {"wavelength_offset_um": 0., "scale": 1., "background": 0., "fwhm_um": 0.}
        for path, value in zip(paths, values):
            parts = path.split(".")
            if len(parts) == 3 and parts[0] == "dataset" and int(parts[1]) == dataset_index:
                result[parts[2]] = float(value)
        return result

    def predict_dataset(current, dataset_index, wavelengths, channels, polarization, values):
        settings = nuisance(values, dataset_index)
        shifted = wavelengths + settings["wavelength_offset_um"]
        predicted = _predict(replace(current, polarization=polarization), shifted, channels)
        predicted = _instrument_response(predicted, wavelengths, settings["fwhm_um"])
        return settings["background"] + settings["scale"]*predicted

    def residual(values, selected="train", synthetic=None):
        current = build(values)
        pieces = []
        for d, (wavelengths, channels, measured, pol, train) in enumerate(parsed):
            truth = measured if synthetic is None else synthetic[d]
            predicted = predict_dataset(current, d, wavelengths, channels, pol, values)
            mask = train if selected == "train" else ~train
            if selected == "all": mask = np.ones(len(wavelengths), dtype=bool)
            if mask.any(): pieces.append((predicted[mask]-truth[mask]).ravel())
        return np.concatenate(pieces)

    x0 = np.asarray(initial)
    if method == "global_then_local":
        global_result = differential_evolution(lambda x: float(np.mean(residual(x)**2)),
            list(zip(lower, upper)), seed=seed, maxiter=18, popsize=7, polish=False,
            updating="immediate")
        x0 = global_result.x
    fit = least_squares(residual, x0, bounds=(lower, upper), max_nfev=160,
                        xtol=1e-10, ftol=1e-10, gtol=1e-10, x_scale="jac")
    train_residual = residual(fit.x, "train")
    validation_residual = residual(fit.x, "validation")
    all_residual = residual(fit.x, "all")
    dof = max(1, len(train_residual)-len(paths))
    covariance = np.linalg.pinv(fit.jac.T@fit.jac) * float(np.sum(train_residual**2)/dof)
    sigma = np.sqrt(np.maximum(0, np.diag(covariance)))
    denom = np.outer(sigma, sigma)
    correlation = np.divide(covariance, denom, out=np.zeros_like(covariance), where=denom>0)
    estimates = [{"path": path, "label": _label(path), "initial": initial[i],
                  "value": float(fit.x[i]), "lower": lower[i], "upper": upper[i],
                  "sigma": float(sigma[i]), "ci95": [float(fit.x[i]-1.96*sigma[i]),
                                                       float(fit.x[i]+1.96*sigma[i])],
                  "at_bound": bool(abs(fit.x[i]-lower[i]) < 1e-5*(upper[i]-lower[i]) or
                                   abs(fit.x[i]-upper[i]) < 1e-5*(upper[i]-lower[i]))}
                 for i, path in enumerate(paths)]
    rows = []
    fitted_model = build(fit.x)
    for d, (wavelengths, channels, measured, pol, train) in enumerate(parsed):
        predicted = predict_dataset(fitted_model, d, wavelengths, channels, pol, fit.x)
        for i, wavelength in enumerate(wavelengths):
            rows.append({"dataset": d+1, "polarization": pol,
                         "wavelength_um": float(wavelength), "set": "train" if train[i] else "validation",
                         "measured": {key: float(measured[i,j]) for j,key in enumerate(channels)},
                         "simulated": {key: float(predicted[i,j]) for j,key in enumerate(channels)},
                         "residual": {key: float(predicted[i,j]-measured[i,j]) for j,key in enumerate(channels)}})
    bootstrap_values = []
    if bootstrap:
        rng = np.random.default_rng(seed)
        base_predictions = []
        base_residuals = []
        for d, (wavelengths, channels, measured, pol, _) in enumerate(parsed):
            pred = predict_dataset(fitted_model, d, wavelengths, channels, pol, fit.x)
            base_predictions.append(pred); base_residuals.append(measured-pred)
        for bootstrap_index in range(bootstrap):
            progress(bootstrap_index, bootstrap, "Residual bootstrap fits")
            synthetic = []
            for prediction, errors in zip(base_predictions, base_residuals):
                choice = rng.integers(0, len(errors), len(errors))
                synthetic.append(np.clip(prediction+errors[choice], 0, 1))
            trial = least_squares(lambda x: residual(x, "train", synthetic), fit.x,
                bounds=(lower, upper), max_nfev=70, x_scale="jac")
            bootstrap_values.append(trial.x.tolist())
    identifiability = []
    for i, estimate in enumerate(estimates):
        relative_ci = 3.92*sigma[i]/max(abs(fit.x[i]), 1e-12)
        max_corr = max([abs(correlation[i,j]) for j in range(len(paths)) if j != i] or [0])
        identifiability.append({"path": estimate["path"], "relative_ci95": float(relative_ci),
            "maximum_absolute_correlation": float(max_corr),
            "status": "weak" if relative_ci > .5 or max_corr > .95 or estimate["at_bound"] else "acceptable"})
    return {"parameters": estimates, "correlation": correlation.tolist(),
            "identifiability": identifiability, "rows": rows,
            "train_rmse": float(np.sqrt(np.mean(train_residual**2))),
            "validation_rmse": None if not len(validation_residual) else float(np.sqrt(np.mean(validation_residual**2))),
            "all_rmse": float(np.sqrt(np.mean(all_residual**2))),
            "success": bool(fit.success), "message": fit.message, "evaluations": int(fit.nfev),
            "method": method, "bootstrap_values": bootstrap_values,
            "nuisance_convention": "Simulation is evaluated at measured wavelength + fitted offset, convolved with a Gaussian instrument FWHM, then transformed as background + scale × signal.",
            "note": "Linearized covariance and residual bootstrap are conditional on the selected model, bounds, noise structure, and optimizer basin."}
