"""Field and power comparison at resonance and reference wavelengths."""

from dataclasses import asdict

import numpy as np

from stack import StackModel, solve_stack, vertical_field


def resonance_field_report(model: StackModel, center_um: float, linewidth_um: float,
                           plane: str = "xz", lateral_points: int = 41,
                           points_per_layer: int = 15) -> dict:
    if not np.isfinite([center_um, linewidth_um]).all() or center_um <= 0 or linewidth_um <= 0:
        raise ValueError("Resonance center and linewidth must be finite and positive")
    if linewidth_um >= center_um/2:
        raise ValueError("Linewidth is too large for a local resonance field comparison")
    wavelengths = [("lower half-width", center_um-linewidth_um/2),
                   ("fitted center", center_um),
                   ("upper half-width", center_um+linewidth_um/2),
                   ("off resonance", center_um+3*linewidth_um)]
    cases = []
    for label, wavelength in wavelengths:
        current = StackModel(**{**asdict(model), "layers": model.layers,
                                "wavelength_um": float(wavelength)})
        field = vertical_field(current, plane, .5, lateral_points, points_per_layer)
        scattering = solve_stack(current)
        cases.append({"label": label, "wavelength_um": float(wavelength),
                      "R": scattering["R"], "T": scattering["T"], "A": scattering["A"],
                      "max_E2": float(np.max(field["E2"])),
                      "mean_E2": float(np.mean(field["E2"])),
                      "layer_absorption": field["layer_absorption"],
                      "field": field})
    center = cases[1]
    off = cases[3]
    return {"cases": cases, "center_to_off_max_E2": center["max_E2"]/max(off["max_E2"], 1e-30),
            "center_to_off_mean_E2": center["mean_E2"]/max(off["mean_E2"], 1e-30),
            "warning": "The comparison uses a fitted far-field linewidth and fixed spatial sampling. Confirm the fit, Fourier/grid convergence, and field-sampling convergence before attributing localization to one mode."}
