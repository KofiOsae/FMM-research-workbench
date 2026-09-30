"""Transparent spectral metrics for ring and knot resonators."""

from __future__ import annotations
import math
import numpy as np

C_UM_S = 299_792_458.0 * 1e6


def resonator_metrics(wavelength_um: float, linewidth_um: float,
                      group_index: float, round_trip_length_um: float,
                      measured_fsr_um: float | None = None,
                      transmission_minimum: float | None = None) -> dict:
    values = [wavelength_um, linewidth_um, group_index, round_trip_length_um]
    if not np.isfinite(values).all() or min(values) <= 0:
        raise ValueError("Wavelength, linewidth, group index, and round-trip length must be positive")
    if linewidth_um >= wavelength_um:
        raise ValueError("Linewidth must be smaller than the center wavelength")
    estimated_fsr = wavelength_um**2 / (group_index * round_trip_length_um)
    fsr = estimated_fsr if measured_fsr_um is None else float(measured_fsr_um)
    if not np.isfinite(fsr) or fsr <= 0:
        raise ValueError("Measured FSR must be positive")
    loaded_q = wavelength_um / linewidth_um
    frequency_hz = C_UM_S / wavelength_um
    omega = 2 * math.pi * frequency_hz
    photon_lifetime_s = loaded_q / omega
    round_trip_time_s = group_index * round_trip_length_um / C_UM_S
    result = {
        "loaded_q": loaded_q, "estimated_fsr_um": estimated_fsr,
        "used_fsr_um": fsr, "finesse": fsr / linewidth_um,
        "photon_lifetime_s": photon_lifetime_s,
        "round_trip_time_s": round_trip_time_s,
        "lifetime_round_trips": photon_lifetime_s / round_trip_time_s,
        "round_trip_length_um": round_trip_length_um,
        "equivalent_ring_diameter_um": round_trip_length_um / math.pi,
        "scope": "These are spectral/circuit metrics. They do not calculate bend radiation, contact-region coupling, polarization mixing, or mechanical deformation in a microfiber knot."
    }
    if transmission_minimum is not None:
        tmin = float(transmission_minimum)
        if not np.isfinite(tmin) or not 0 <= tmin <= 1:
            raise ValueError("On-resonance transmission must be between 0 and 1")
        root = math.sqrt(tmin)
        total_rate = 1.0 / loaded_q
        candidates = []
        for label, sign in (("intrinsic-loss-dominant branch", 1),
                            ("coupling-dominant branch", -1)):
            intrinsic_rate = total_rate * (1 + sign * root) / 2
            coupling_rate = total_rate * (1 - sign * root) / 2
            candidates.append({"branch": label,
                               "intrinsic_q": math.inf if intrinsic_rate == 0 else 1/intrinsic_rate,
                               "coupling_q": math.inf if coupling_rate == 0 else 1/coupling_rate})
        result["coupling_candidates"] = candidates
        result["coupling_note"] = "Intensity depth alone cannot distinguish the under-coupled and over-coupled branches; phase or an independent loss measurement is required."
    return result
