"""Fabrication-tolerance propagation through finite-stack spectra."""

from dataclasses import asdict

import numpy as np

from stack import StackModel, solve_stack
from sweep import set_parameter


def tolerance_study(model: StackModel, parameters: list[dict], samples: int,
                    wavelength_start: float, wavelength_stop: float,
                    wavelength_points: int, quantity: str = "R",
                    extremum: str = "max", seed: int = 12345,
                    correlation: list[list[float]] | None = None) -> dict:
    if not 5 <= samples <= 200:
        raise ValueError("Use 5–200 tolerance samples")
    if not 7 <= wavelength_points <= 101:
        raise ValueError("Use 7–101 wavelengths per sample")
    if not 0 < wavelength_start < wavelength_stop or not np.isfinite([wavelength_start, wavelength_stop]).all():
        raise ValueError("Use positive increasing wavelength limits")
    if quantity not in ("R", "T", "A", "R0", "T0") or extremum not in ("max", "min"):
        raise ValueError("Choose a supported quantity and maximum or minimum")
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
    rows, failures = [], []
    for sample in range(samples):
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
                    spectrum.append((float(wavelength), solve_stack(point)[quantity]))
                except (ValueError, np.linalg.LinAlgError):
                    continue
            if len(spectrum) < 3:
                raise ValueError("Fewer than three wavelength points completed")
            values_array = np.asarray([value for _, value in spectrum])
            index = int(np.argmax(values_array) if extremum == "max" else np.argmin(values_array))
            rows.append({"sample": sample+1, "parameters": values,
                         "feature_wavelength_um": spectrum[index][0],
                         "feature_value": float(spectrum[index][1]),
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
    sensitivities = []
    for item in parameters:
        x = np.asarray([row["parameters"][item["path"]] for row in rows])
        for target, y in (("feature_wavelength_um", feature_wavelength),
                          ("feature_value", feature_value)):
            slope = float(np.polyfit(x, y, 1)[0]) if np.std(x) > 0 else 0.0
            correlation = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else 0.0
            sensitivities.append({"parameter": item["path"], "target": target,
                                  "linear_slope": slope, "pearson_r": correlation})
    return {"rows": rows, "failures": failures, "summary": summary,
            "sensitivities": sensitivities, "quantity": quantity, "extremum": extremum,
            "seed": seed, "wavelength_points": wavelength_points,
            "input_correlation": corr.tolist(),
            "empirical_parameter_correlation": np.corrcoef(np.asarray([[row["parameters"][item["path"]] for item in parameters] for row in rows]), rowvar=False).tolist() if len(parameters) > 1 else [[1.0]],
            "warning": "Feature locations are selected from the entered uniform wavelength grid. Refine the grid or use adaptive fitting before interpreting linewidth-scale shifts; distributions describe the entered parameter model only."}
