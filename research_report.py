"""Reproducible Markdown methods report for a finite-stack calculation."""

from datetime import datetime, timezone

from materials import MATERIAL_INFO
from stack import StackModel


def research_report(model: StackModel, software: dict, title: str = "Photonics simulation") -> str:
    model.validate()
    lines = [f"# {title.strip() or 'Photonics simulation'}", "",
             f"Generated: {datetime.now(timezone.utc).isoformat()}", "",
             "## Physical model", "",
             f"- Vacuum wavelength: {model.wavelength_um:.9g} µm",
             f"- Incidence: θ={model.theta_deg:.9g}°, φ={model.phi_deg:.9g}°",
             f"- Polarization: {model.polarization}",
             f"- Incident / exit refractive index: {model.incident_n:.9g} / {model.exit_n:.9g}",
             f"- Primitive periods: {model.period_x_um:.9g} µm, {model.period_y_um:.9g} µm",
             f"- Primitive-vector angle: {model.lattice_angle_deg:.9g}°", "",
             "## Finite layers, physical incident-to-exit order", ""]
    for index, layer in enumerate(model.layers, 1):
        lines.append(f"### Layer {index}")
        lines.extend([f"- Geometry: {layer.kind}", f"- Thickness: {layer.thickness_um:.9g} µm",
                      f"- Background: {layer.background_material}; constant n={layer.background_n:.9g}",
                      f"- Feature: {layer.feature_material}; constant n={layer.feature_n:.9g}",
                      f"- Fill parameters: x={layer.fill_x:.9g}, y={layer.fill_y:.9g}",
                      f"- Feature offsets: a1={layer.offset_x:.9g}, a2={layer.offset_y:.9g}", ""])
    keys = sorted({value for layer in model.layers
                   for value in (layer.background_material, layer.feature_material)
                   if value not in ("dielectric", "air")})
    lines.extend(["## Material provenance", ""])
    if not keys:
        lines.append("- Air and/or explicitly constant refractive indices were used.")
    for key in keys:
        item = MATERIAL_INFO.get(key)
        if item:
            lines.append(f"- {item[0]}: {item[2]}; range {item[1]}.")
        else:
            lines.append(f"- {key}: see the locally imported material record and retain its source file.")
    lines.extend(["", "## Numerical settings", "",
                  f"- Requested FMM harmonics: {model.order_budget}",
                  f"- Geometry grid: {model.grid_size} samples per primitive axis",
                  f"- Software: " + ", ".join(f"{key} {value}" for key, value in software.items()), "",
                  "## Required result-specific record", "",
                  "Add the exported result JSON and record:", "",
                  "- Actual retained Fourier orders and order-to-order changes.",
                  "- Geometry-grid refinement change.",
                  "- R, T, A and the energy-accounting residual.",
                  "- Wavelength/angle sampling and any interpolation used only for display.",
                  "- Fit model, residuals, confidence assumptions, validation data, and parameter bounds.",
                  "- Field sampling, reference plane, component, normalization, and layer-loss mismatch.", "",
                  "## Interpretation limits", "",
                  "- Numerical convergence does not include fabrication or material uncertainty.",
                  "- A far-field spectral feature alone does not prove a guided mode, BIC, or band gap.",
                  "- Homogenized slab phase matching establishes momentum compatibility only.",
                  "- Driven Jones/Stokes maps include the nonresonant scattering background.",
                  "- Compare the final claim with an independent solver or measurement.", ""])
    return "\n".join(lines)
