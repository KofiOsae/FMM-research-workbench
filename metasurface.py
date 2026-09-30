"""Local-periodic metasurface phase library and radial flat-lens assignment."""

from __future__ import annotations

from dataclasses import asdict
import numpy as np

import run_jobs
from stack import StackModel, solve_stack
from sweep import set_parameter, parameter_label, sweep_values


def phase_library_and_lens(model: StackModel, parameter: str, start: float,
                           stop: float, points: int, port: str,
                           lens_radius_um: float, focal_length_um: float,
                           radial_points: int, minimum_power: float = 0.0) -> dict:
    """Build a specular phase library, then assign it to a radial lens profile.

    This is a local periodic approximation (LPA): every library point is an
    infinite periodic array. It does not propagate the assembled finite lens.
    """
    if model.polarization not in ("s", "p"):
        raise ValueError("A coherent s or p polarization is required for a phase library")
    if abs(model.theta_deg) > 1e-9:
        raise ValueError("The flat-lens assistant currently requires normal incidence")
    if port not in ("transmitted", "reflected"):
        raise ValueError("Phase-library port must be transmitted or reflected")
    if not np.isfinite([lens_radius_um, focal_length_um, minimum_power]).all() \
            or lens_radius_um <= 0 or focal_length_um <= 0 or not 0 <= minimum_power <= 1:
        raise ValueError("Lens radius and focal length must be positive; minimum power is 0–1")
    if not 3 <= radial_points <= 2001:
        raise ValueError("Use 3–2001 radial lens samples")
    values = sweep_values(start, stop, points)
    phase_key = "t_phase_deg" if port == "transmitted" else "r_phase_deg"
    power_key = "T0" if port == "transmitted" else "R0"
    rows = []
    for index, value in enumerate(values):
        run_jobs.progress(index, len(values) + radial_points, "Metasurface phase library")
        current = set_parameter(model, parameter, float(value))
        result = solve_stack(current)
        rows.append({"parameter": float(value), "phase_deg": result[phase_key],
                     "power": float(result[power_key]), "R": result["R"],
                     "T": result["T"], "A": result["A"],
                     "actual_orders": result["actual_orders"]})
    usable = [row for row in rows if row["phase_deg"] is not None and row["power"] >= minimum_power]
    if len(usable) < 2:
        raise ValueError("Fewer than two phase-library entries meet the power threshold")

    radius = np.linspace(0.0, lens_radius_um, radial_points)
    target = (-360.0 / model.wavelength_um *
              (np.sqrt(focal_length_um**2 + radius**2) - focal_length_um)) % 360.0
    phases = np.asarray([row["phase_deg"] for row in usable])
    # A global phase does not change the focus. Search it to minimize the
    # quantization error produced by this finite library.
    best = None
    for offset in np.linspace(0, 360, 721, endpoint=False):
        errors = np.angle(np.exp(1j*np.deg2rad(phases[:, None] - target[None, :] - offset)), deg=True)
        choices = np.argmin(np.abs(errors), axis=0)
        selected = errors[choices, np.arange(radial_points)]
        score = float(np.sqrt(np.mean(selected**2)))
        if best is None or score < best[0]:
            best = (score, float(offset), choices, selected)
    rms, offset, choices, errors = best
    layout = []
    for index, (r, wanted, choice, error) in enumerate(zip(radius, target, choices, errors)):
        run_jobs.progress(len(values) + index, len(values) + radial_points, "Assigning flat-lens cells")
        chosen = usable[int(choice)]
        layout.append({"radius_um": float(r), "target_phase_deg": float(wanted),
                       "global_offset_deg": offset, "parameter": chosen["parameter"],
                       "library_phase_deg": chosen["phase_deg"],
                       "phase_error_deg": float(error), "power": chosen["power"]})
    sorted_phase = np.sort(np.asarray([row["phase_deg"] % 360 for row in usable]))
    gaps = np.diff(np.r_[sorted_phase, sorted_phase[0] + 360])
    numerical_aperture = model.exit_n * lens_radius_um / np.sqrt(
        lens_radius_um**2 + focal_length_um**2)
    return {
        "library": rows, "layout": layout, "parameter": parameter,
        "parameter_label": parameter_label(parameter), "port": port,
        "wavelength_um": model.wavelength_um, "lens_radius_um": lens_radius_um,
        "focal_length_um": focal_length_um, "numerical_aperture": float(numerical_aperture),
        "rms_phase_error_deg": rms, "maximum_phase_gap_deg": float(np.max(gaps)),
        "mean_selected_power": float(np.mean([row["power"] for row in layout])),
        "minimum_selected_power": float(np.min([row["power"] for row in layout])),
        "global_phase_offset_deg": offset, "minimum_library_power": minimum_power,
        "model": asdict(model),
        "scope": "Local periodic approximation: library cells are infinite periodic arrays. The radial assignment is not a finite-aperture propagation solve.",
        "validation": "Repeat the library with larger Fourier budget and finer geometry grid; validate the assembled device with a finite full-wave or propagation model, especially for high NA or rapidly varying cells."
    }
