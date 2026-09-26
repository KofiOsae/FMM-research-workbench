"""Approximate resonant Jones/Stokes separation from complex spectral sidebands."""

from __future__ import annotations

from dataclasses import replace
import numpy as np

from stack import StackModel, solve_stack


def _complex(item):
    return complex(float(item["real"]), float(item["imag"]))


def _stokes(ep: complex, es: complex) -> dict:
    s0 = float(abs(ep)**2+abs(es)**2)
    s1 = float(abs(ep)**2-abs(es)**2)
    s2 = float(2*np.real(ep*np.conj(es)))
    s3 = float(-2*np.imag(ep*np.conj(es)))
    return {"Ep": {"real": float(ep.real), "imag": float(ep.imag)},
            "Es": {"real": float(es.real), "imag": float(es.imag)},
            "S0": s0, "S1": s1, "S2": s2, "S3": s3,
            "normalized": {"S1": s1/max(s0, 1e-30),
                           "S2": s2/max(s0, 1e-30),
                           "S3": s3/max(s0, 1e-30)},
            "orientation_deg": float(.5*np.degrees(np.arctan2(s2, s1))) if s0 else 0.,
            "ellipticity_deg": float(.5*np.degrees(np.arcsin(np.clip(s3/max(s0, 1e-30), -1, 1))))}


def _channel(model, wavelength, m, n, port):
    result = solve_stack(replace(model, wavelength_um=float(wavelength)))
    row = next((item for item in result.get("order_amplitudes", [])
                if item["m"] == m and item["n"] == n), None)
    if row is None:
        raise ValueError(f"Diffraction order ({m},{n}) is not propagating at {wavelength:.7g} µm")
    item = row[f"{port}_polarization"]
    return _complex(item["Ep"]), _complex(item["Es"])


def resonant_polarization(model: StackModel, center_um: float, linewidth_um: float,
                          m: int = 0, n: int = 0, port: str = "reflected",
                          sideband_span: float = 3.) -> dict:
    """Subtract a locally linear complex background around a fitted resonance."""
    model.validate()
    if model.polarization == "unpolarized":
        raise ValueError("Resonant Jones analysis requires coherent s or p input")
    if port not in ("reflected", "transmitted") or linewidth_um <= 0 or sideband_span < 1.5:
        raise ValueError("Use reflected/transmitted, positive linewidth, and sideband span ≥ 1.5")
    if center_um-sideband_span*linewidth_um <= 0:
        raise ValueError("Sideband wavelength must remain positive")

    center = _channel(model, center_um, m, n, port)
    estimates = []
    for span in (sideband_span, sideband_span+1):
        low = _channel(model, center_um-span*linewidth_um, m, n, port)
        high = _channel(model, center_um+span*linewidth_um, m, n, port)
        background = tuple((a+b)/2 for a, b in zip(low, high))
        resonant = tuple(a-b for a, b in zip(center, background))
        estimates.append({"span_linewidths": float(span),
                          "low_wavelength_um": float(center_um-span*linewidth_um),
                          "high_wavelength_um": float(center_um+span*linewidth_um),
                          "background": _stokes(*background),
                          "resonant_residual": _stokes(*resonant)})
    a, b = estimates[0]["resonant_residual"], estimates[1]["resonant_residual"]
    ja = np.asarray([complex(**a[key]) for key in ("Ep", "Es")])
    jb = np.asarray([complex(**b[key]) for key in ("Ep", "Es")])
    relative_change = float(np.linalg.norm(ja-jb)/max(np.linalg.norm(ja), 1e-30))
    return {"center_um": float(center_um), "linewidth_um": float(linewidth_um),
            "order": {"m": int(m), "n": int(n)}, "port": port,
            "driven_center": _stokes(*center), "estimates": estimates,
            "sideband_stability_relative_change": relative_change,
            "interpretation": "A locally linear complex background was interpolated from symmetric spectral sidebands and subtracted from the driven center field.",
            "limitation": "This is a sideband decomposition, not an S-matrix pole residue or eigenmode. Report the sideband span and stability; do not use it alone to claim a BIC charge."}
