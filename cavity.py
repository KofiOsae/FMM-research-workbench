"""Single-mode cavity Purcell estimate; not a dipole-field simulation."""

import math


def purcell_estimate(cavity_wavelength_um: float, emitter_wavelength_um: float,
                     refractive_index: float, quality_factor: float,
                     mode_volume_um3: float, overlap: float) -> dict:
    values = (cavity_wavelength_um, emitter_wavelength_um, refractive_index,
              quality_factor, mode_volume_um3, overlap)
    if not all(math.isfinite(v) for v in values) or min(values[:5]) <= 0 or not 0 <= overlap <= 1:
        raise ValueError("Use positive finite wavelengths, n, Q, V and overlap between 0 and 1")
    normalized_volume = mode_volume_um3/(cavity_wavelength_um/refractive_index)**3
    maximum = 3/(4*math.pi**2)*quality_factor/normalized_volume
    relative_frequency = cavity_wavelength_um/emitter_wavelength_um
    spectral_factor = 1/(1+4*quality_factor**2*(relative_frequency-1)**2)
    return {"normalized_mode_volume": normalized_volume,
            "ideal_on_resonance": maximum,
            "spectral_factor": spectral_factor,
            "estimated_mode_contribution": maximum*spectral_factor*overlap,
            "method": "single-mode weak-coupling Q/V estimate",
            "note": "Requires independently established Q and electromagnetic mode volume; overlap combines position and dipole orientation. This is not a total dipole decay-rate or a structure-specific FMM result."}
