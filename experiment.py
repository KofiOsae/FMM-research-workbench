"""Measured spectrum comparison and bounded uniform-film thickness fit."""

from dataclasses import asdict
import csv
import io

import numpy as np
from scipy.optimize import least_squares

from stack import StackModel, Layer
from tmm import solve_tmm


def parse_measurement(csv_text: str) -> tuple[np.ndarray, list[str], np.ndarray]:
    if not isinstance(csv_text, str) or len(csv_text) > 100_000:
        raise ValueError("Measurement CSV is empty or exceeds 100 kB")
    reader = csv.DictReader(io.StringIO(csv_text))
    header = reader.fieldnames or []
    channels = [key for key in ("R", "T", "A") if key in header]
    if "wavelength_um" not in header or not channels or any(key not in ("wavelength_um", "R", "T", "A") for key in header):
        raise ValueError("CSV needs wavelength_um and at least one of R,T,A; power values are fractions 0–1")
    rows = []
    for row in reader:
        if len(rows) >= 200:
            raise ValueError("Use at most 200 measured wavelengths")
        try:
            values = [float(row[key]) for key in ["wavelength_um", *channels]]
        except (ValueError, TypeError, KeyError) as exc:
            raise ValueError("Every measured row must have numeric values") from exc
        if not np.isfinite(values).all() or values[0] <= 0 or any(not 0 <= v <= 1 for v in values[1:]):
            raise ValueError("Require finite wavelength > 0 and measured R/T/A fractions from 0 to 1")
        if rows and values[0] <= rows[-1][0]:
            raise ValueError("Measured wavelengths must increase strictly")
        rows.append(values)
    if len(rows) < 5:
        raise ValueError("Provide at least five measured wavelengths")
    array = np.asarray(rows)
    if all(key in channels for key in ("R", "T", "A")):
        sums = array[:, 1:].sum(axis=1)
        if np.max(np.abs(sums-1)) > .05:
            raise ValueError("Measured R+T+A differs from one by more than 0.05")
    return array[:,0], channels, array[:,1:]


def compare_or_fit(model: StackModel, csv_text: str, layer_index: int,
                   fit: bool, lower_um: float = .01, upper_um: float = 2.0) -> dict:
    model.validate()
    if any(layer.kind != "uniform" for layer in model.layers):
        raise ValueError("Measured-data thickness fitting currently requires uniform layers (transfer matrix)")
    wavelengths, channels, measured = parse_measurement(csv_text)
    if not 0 <= layer_index < len(model.layers):
        raise ValueError("Choose an existing layer to fit")
    if not np.isfinite([lower_um, upper_um]).all() or not 0 < lower_um < upper_um:
        raise ValueError("Thickness bounds must be positive and increasing")

    def prediction(thickness: float) -> np.ndarray:
        layers = list(model.layers)
        layers[layer_index] = Layer(**{**asdict(layers[layer_index]), "thickness_um": thickness})
        rows = []
        for wavelength in wavelengths:
            current = StackModel(**{**asdict(model), "layers": tuple(layers),
                                    "wavelength_um": float(wavelength)})
            result = solve_tmm(current)
            rows.append([result[key] for key in channels])
        return np.asarray(rows)

    initial = model.layers[layer_index].thickness_um
    if fit:
        if not lower_um <= initial <= upper_um:
            raise ValueError("Initial thickness must lie inside the fitting bounds")
        fit_result = least_squares(lambda value: (prediction(float(value[0]))-measured).ravel(),
                                   [initial], bounds=([lower_um], [upper_um]),
                                   max_nfev=80, xtol=1e-10, ftol=1e-10, gtol=1e-10)
        thickness = float(fit_result.x[0])
    else:
        thickness = initial
    simulated = prediction(thickness)
    residual = simulated-measured
    rows = [{"wavelength_um": float(wavelengths[i]),
             "measured": {key: float(measured[i,j]) for j,key in enumerate(channels)},
             "simulated": {key: float(simulated[i,j]) for j,key in enumerate(channels)},
             "residual": {key: float(residual[i,j]) for j,key in enumerate(channels)}}
            for i in range(len(wavelengths))]
    return {"rows": rows, "channels": channels, "fit": fit,
            "layer_index": layer_index, "thickness_um": thickness,
            "initial_thickness_um": initial, "bounds_um": [lower_um, upper_um],
            "rmse": float(np.sqrt(np.mean(residual**2))),
            "max_absolute_residual": float(np.max(np.abs(residual))),
            "method": "bounded nonlinear least squares; unweighted power-fraction residuals; uniform-stack TMM",
            "note": "A fit is conditional on this model and does not establish a unique physical structure or parameter uncertainty."}
