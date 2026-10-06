"""Fabrication-tolerance propagation through finite-stack spectra."""

from dataclasses import asdict

from run_jobs import progress

import numpy as np

from stack import StackModel, solve_stack
from sweep import set_parameter


def tolerance_study(model: StackModel, parameters: list[dict], samples: int,
                    wavelength_start: float, wavelength_stop: float,
                    wavelength_points: int, quantity: str = "R",
                    extremum: str = "max", seed: int = 12345,
                    correlation: list[list[float]] | None = None,
                    order_m: int = 0, order_n: int = 0,
                    operating_wavelength_um: float | None = None) -> dict:
    if not 5 <= samples <= 200:
        raise ValueError("Use 5–200 tolerance samples")
    if not 7 <= wavelength_points <= 101:
        raise ValueError("Use 7–101 wavelengths per sample")
    if not 0 < wavelength_start < wavelength_stop or not np.isfinite([wavelength_start, wavelength_stop]).all():
        raise ValueError("Use positive increasing wavelength limits")
    if quantity not in ("R", "T", "A", "R0", "T0", "R_order", "T_order") or extremum not in ("max", "min"):
        raise ValueError("Choose a supported quantity and maximum or minimum")
    if operating_wavelength_um is not None and (not np.isfinite(operating_wavelength_um) or operating_wavelength_um <= 0):
        raise ValueError("Operating wavelength must be positive and finite")
    if not 1 <= len(parameters) <= 5:
        raise ValueError("Use 1–5 uncertain parameters")
    for item in parameters:
        if not item.get("path") or not np.isfinite(float(item.get("sigma", 0))) or float(item["sigma"]) <= 0:
            raise ValueError("Every uncertain parameter needs a path and positive sigma")
    rng = np.random.default_rng(seed)
    count = len(parameters)
    corr = np.eye(count) if correlation is None else np.asarray(correlation, dtype=float)
    if corr.shape != (count, count) or not np.isfinite(corr).all():
        raise ValueError(f"Correlation matrix must be {count} by {count} and finite")
    if not np.allclose(corr, corr.T, atol=1e-10) or not np.allclose(np.diag(corr), 1, atol=1e-10):
        raise ValueError("Correlation matrix must be symmetric with ones on its diagonal")
    eigenvalues = np.linalg.eigvalsh(corr)
    if eigenvalues[0] < -1e-10:
        raise ValueError("Correlation matrix must be positive semidefinite")
    # Eigen construction also supports exactly correlated parameters, where a
    # Cholesky factor would fail.
    factor = np.linalg.eigh(corr)
    transform = factor[1] @ np.diag(np.sqrt(np.maximum(factor[0], 0)))
    standard_draws = rng.normal(size=(samples, count)) @ transform.T
    wavelengths = np.linspace(wavelength_start, wavelength_stop, wavelength_points)
    def response(result: dict) -> float:
        if quantity in ("R", "T", "A", "R0", "T0"):
            return float(result[quantity])
        port = quantity[0]
        item = next((row for row in result["orders"]
                     if row["m"] == order_m and row["n"] == order_n), None)
        return 0.0 if item is None else float(item[port])
    rows, failures = [], []
    for sample in range(samples):
        progress(sample, samples, "Fabrication tolerance samples")
        current = model
        values = {}
        try:
            for parameter_index, item in enumerate(parameters):
                value = float(item["mean"]) + float(item["sigma"])*standard_draws[sample, parameter_index]
                if "lower" in item:
                    value = max(value, float(item["lower"]))
                if "upper" in item:
                    value = min(value, float(item["upper"]))
                current = set_parameter(current, str(item["path"]), value)
                values[item["path"]] = value
            spectrum = []
            for wavelength in wavelengths:
                point = StackModel(**{**asdict(current), "layers": current.layers,
                                      "wavelength_um": float(wavelength)})
                try:
                    spectrum.append((float(wavelength), response(solve_stack(point))))
                except (ValueError, np.linalg.LinAlgError):
                    continue
            if len(spectrum) < 3:
                raise ValueError("Fewer than three wavelength points completed")
            values_array = np.asarray([value for _, value in spectrum])
            index = int(np.argmax(values_array) if extremum == "max" else np.argmin(values_array))
            operating_value = None
            if operating_wavelength_um is not None:
                operating_model = StackModel(**{**asdict(current), "layers": current.layers,
                                                "wavelength_um": float(operating_wavelength_um)})
                operating_value = response(solve_stack(operating_model))
            rows.append({"sample": sample+1, "parameters": values,
                         "feature_wavelength_um": spectrum[index][0],
                         "feature_value": float(spectrum[index][1]),
                         "operating_value": operating_value,
                         "calculated_wavelength_points": len(spectrum)})
        except (ValueError, np.linalg.LinAlgError) as exc:
            failures.append({"sample": sample+1, "parameters": values, "reason": str(exc)})
    if len(rows) < 3:
        reason = failures[0]["reason"] if failures else "unknown numerical failure"
        raise ValueError(f"Fewer than three tolerance samples completed; first failure: {reason}")
    feature_wavelength = np.asarray([row["feature_wavelength_um"] for row in rows])
    feature_value = np.asarray([row["feature_value"] for row in rows])
    summary = {}
    for name, values in (("feature_wavelength_um", feature_wavelength),
                         ("feature_value", feature_value)):
        summary[name] = {"mean": float(np.mean(values)), "std": float(np.std(values, ddof=1)),
                         "p05": float(np.percentile(values, 5)),
                         "p50": float(np.percentile(values, 50)),
                         "p95": float(np.percentile(values, 95))}
    operating_values = None
    if operating_wavelength_um is not None:
        operating_values = np.asarray([row["operating_value"] for row in rows])
        summary["operating_value"] = {"mean": float(np.mean(operating_values)),
                                      "std": float(np.std(operating_values, ddof=1)),
                                      "p05": float(np.percentile(operating_values, 5)),
                                      "p50": float(np.percentile(operating_values, 50)),
                                      "p95": float(np.percentile(operating_values, 95))}
    sensitivities = []
    sensitivity_targets = [("feature_wavelength_um", feature_wavelength),
                           ("feature_value", feature_value)]
    if operating_values is not None:
        sensitivity_targets.append(("operating_value", operating_values))
    for item in parameters:
        x = np.asarray([row["parameters"][item["path"]] for row in rows])
        for target, y in sensitivity_targets:
            slope = float(np.polyfit(x, y, 1)[0]) if np.std(x) > 0 else 0.0
            correlation = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else 0.0
            sensitivities.append({"parameter": item["path"], "target": target,
                                  "linear_slope": slope, "pearson_r": correlation})
    return {"rows": rows, "failures": failures, "summary": summary,
            "sensitivities": sensitivities, "quantity": quantity, "extremum": extremum,
            "order_m": int(order_m), "order_n": int(order_n),
            "operating_wavelength_um": operating_wavelength_um,
            "seed": seed, "wavelength_points": wavelength_points,
            "input_correlation": corr.tolist(),
            "empirical_parameter_correlation": np.corrcoef(np.asarray([[row["parameters"][item["path"]] for item in parameters] for row in rows]), rowvar=False).tolist() if len(parameters) > 1 else [[1.0]],
            "warning": "Feature locations are selected from the entered uniform wavelength grid. A missing diffraction order has zero far-field efficiency at that point. Refine the grid before interpreting wavelength-scale shifts; distributions describe the entered parameter model only."}
