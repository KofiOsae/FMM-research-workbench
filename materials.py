"""Measured optical constants with explicit wavelength coverage (µm)."""

import csv
import json
import re
import threading
import os
from functools import lru_cache
from pathlib import Path

import numpy as np


DATA = Path(__file__).with_name("gold_johnson_christy.csv")
SILICON_DATA = Path(__file__).with_name("silicon_green_2008.yml")
USER_DATA = Path(__file__).with_name("user_materials.json")
_material_lock = threading.RLock()
MATERIAL_INFO = {
    "air": ("Air (n = 1 approximation)", "All model wavelengths", "Idealized constant index"),
    "silica": ("Fused silica", "0.21–6.7 µm", "Malitson 1965 Sellmeier; 20 °C"),
    "bk7": ("N-BK7 optical glass", "0.3–2.5 µm", "SCHOTT 2017 Sellmeier; 20 °C"),
    "n_sf11": ("N-SF11 dense flint optical glass", "0.365–2.5 µm", "SCHOTT optical-glass data sheet; coefficients also indexed by refractiveindex.info"),
    "n_bak4": ("N-BAK4 barium crown optical glass", "0.365–2.5 µm", "SCHOTT optical-glass data sheet; coefficients also indexed by refractiveindex.info"),
    "caf2": ("Calcium fluoride", "0.23–9.7 µm", "Malitson 1963 Sellmeier; 24 °C"),
    "silicon": ("Crystalline silicon", "0.25–1.45 µm", "Green 2008 n,k; 300 K"),
    "gold": ("Gold", "0.1879–1.9370 µm", "Johnson and Christy 1972 n,k"),
    "baf2": ("Barium fluoride", "0.2652–10.346 µm", "Malitson 1964 Sellmeier; 25 °C, k omitted"),
    "lif": ("Lithium fluoride", "0.10–11 µm", "Li 1976 Sellmeier; 24 °C, k omitted"),
    "kf": ("Potassium fluoride", "0.15–22 µm", "Li 1976 Sellmeier; 24 °C, k omitted"),
    "nacl": ("Sodium chloride", "0.20–30 µm", "Li 1976 Sellmeier; 24 °C, k omitted"),
    "tio2": ("Rutile TiO₂ · ordinary ray", "0.43–1.53 µm", "DeVore 1951; room temperature, k omitted"),
    "al2o3": ("α-Al₂O₃ sapphire · ordinary ray", "0.20–5.0 µm", "Malitson and Dodge 1972; 20 °C, k omitted"),
    "mgf2": ("MgF₂ · ordinary ray", "0.20–7.0 µm", "Dodge 1984; 19 °C, k omitted"),
    "sin": ("Silicon nitride · measured waveguide film", "0.45–1.65 µm", "Sellmeier fit, APL Photonics 2026, DOI 10.1063/5.0293497"),
    "silver_rakic": ("Silver · Lorentz–Drude", "0.207–12.398 µm", "Rakic et al. 1998, DOI 10.1364/AO.37.005271"),
    "aluminum_rakic": ("Aluminum · Lorentz–Drude", "0.207–12.398 µm", "Rakic et al. 1998, DOI 10.1364/AO.37.005271"),
    "copper_rakic": ("Copper · Lorentz–Drude", "0.207–12.398 µm", "Rakic et al. 1998, DOI 10.1364/AO.37.005271"),
    "pmma": ("PMMA · measured bulk polymer", "0.4368–1.052 µm", "Sultanova et al. 2009, DOI 10.12693/APhysPolA.116.585"),
    "polycarbonate": ("Polycarbonate · measured bulk polymer", "0.4368–1.052 µm", "Sultanova et al. 2009, DOI 10.12693/APhysPolA.116.585"),
    "polystyrene": ("Polystyrene · measured bulk polymer", "0.4368–1.052 µm", "Sultanova et al. 2009, DOI 10.12693/APhysPolA.116.585"),
    "zeonex": ("Zeonex E48R · measured bulk polymer", "0.4368–1.052 µm", "Sultanova et al. 2009, DOI 10.12693/APhysPolA.116.585"),
}

POLYMER_WAVELENGTHS = np.asarray([.4368, .4861, .5876, .6328, .703, .833, .879, 1.052])
POLYMER_INDICES = {
    "pmma": [1.502, 1.497, 1.491, 1.489, 1.486, 1.484, 1.483, 1.481],
    "polycarbonate": [1.612, 1.599, 1.585, 1.580, 1.575, 1.569, 1.568, 1.565],
    "polystyrene": [1.617, 1.606, 1.592, 1.587, 1.582, 1.577, 1.576, 1.572],
    "zeonex": [1.543, 1.538, 1.531, 1.528, 1.526, 1.523, 1.522, 1.520],
}

RAKIC_LD = {
    # plasma energy, Drude (strength, damping), Lorentz (strength, damping, energy)
    "silver_rakic": (9.01, (.845, .048), ((.065, 3.886, .816), (.124, .452, 4.481),
        (.011, .065, 8.185), (.840, .916, 9.083), (5.646, 2.419, 20.29))),
    "aluminum_rakic": (14.98, (.523, .047), ((.227, .333, .162), (.050, .312, 1.544),
        (.166, 1.351, 1.808), (.030, 3.382, 3.473))),
    "copper_rakic": (10.83, (.575, .030), ((.061, .378, .291), (.104, 1.056, 2.957),
        (.723, 3.213, 5.300), (.638, 4.305, 11.18))),
}


def _rakic_n(name: str, wavelength_um: float) -> complex:
    energy = 1.239841984/wavelength_um
    if not .1 <= energy <= 6:
        raise ValueError("Rakic Lorentz-Drude models cover photon energies 0.1–6 eV")
    plasma, (f0, gamma0), oscillators = RAKIC_LD[name]
    epsilon = 1-f0*plasma**2/(energy*(energy+1j*gamma0))
    for strength, damping, resonance in oscillators:
        epsilon += strength*plasma**2/(resonance**2-energy**2-1j*damping*energy)
    value = np.sqrt(epsilon)
    return complex(-value if value.imag < 0 else value)


def _sellmeier(wavelength_um: float, lower: float, upper: float,
                offset: float, terms: tuple[tuple[float, float], ...], label: str) -> complex:
    if not np.isfinite(wavelength_um) or not lower <= wavelength_um <= upper:
        raise ValueError(f"{label} formula covers {lower}–{upper} µm")
    square = wavelength_um * wavelength_um
    value = 1 + offset + sum(b*square/(square-c*c) for b, c in terms)
    if value <= 0 or not np.isfinite(value):
        raise ValueError(f"{label} dispersion is invalid at this wavelength")
    return complex(np.sqrt(value))


def _user_materials() -> dict:
    if not USER_DATA.exists():
        return {}
    with _material_lock:
        return json.loads(USER_DATA.read_text(encoding="utf-8"))


def import_material(name: str, source: str, csv_text: str) -> dict:
    """Store a local, explicitly unit-labelled wavelength,n,k table."""
    if not isinstance(name, str) or not 1 <= len(name.strip()) <= 60:
        raise ValueError("Material name must contain 1–60 characters")
    if not isinstance(source, str) or not 1 <= len(source.strip()) <= 300:
        raise ValueError("Give a source, specimen, or measurement note (1–300 characters)")
    if not isinstance(csv_text, str) or len(csv_text) > 60_000:
        raise ValueError("Material CSV must be at most 60 kB")
    slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
    if not slug:
        raise ValueError("Material name must contain letters or digits")
    reader = csv.DictReader(csv_text.splitlines())
    if reader.fieldnames != ["wavelength_um", "n", "k"]:
        raise ValueError("CSV header must be exactly: wavelength_um,n,k")
    values = []
    for row in reader:
        if len(values) >= 5000:
            raise ValueError("Use at most 5000 material data rows")
        try:
            wavelength, n, k = (float(row[key]) for key in ("wavelength_um", "n", "k"))
        except (ValueError, TypeError, KeyError) as exc:
            raise ValueError("Every row needs numeric wavelength_um,n,k") from exc
        if not np.isfinite([wavelength, n, k]).all() or wavelength <= 0 or n <= 0 or k < 0:
            raise ValueError("Require finite wavelength > 0, n > 0, and k ≥ 0")
        if values and wavelength <= values[-1][0]:
            raise ValueError("Wavelengths must increase strictly without duplicates")
        values.append([wavelength, n, k])
    if len(values) < 2:
        raise ValueError("Provide at least two wavelength rows")
    key = "custom:" + slug
    with _material_lock:
        records = _user_materials()
        records[key] = {"name": name.strip(), "source": source.strip(), "rows": values}
        temporary = USER_DATA.with_suffix(".tmp")
        temporary.write_text(json.dumps(records, ensure_ascii=False), encoding="utf-8")
        os.replace(temporary, USER_DATA)
    return {"key": key, "name": name.strip(), "rows": len(values),
            "range": f"{values[0][0]}–{values[-1][0]} µm", "source": source.strip()}


@lru_cache(maxsize=1)
def gold_table() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with DATA.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    return tuple(np.array([float(row[key]) for row in rows])
                 for key in ("wavelength_um", "n", "k"))


def gold_n(wavelength_um: float) -> complex:
    """Interpolate n and k separately; do not extrapolate measured values."""
    wavelengths, n, k = gold_table()
    if not wavelengths[0] <= wavelength_um <= wavelengths[-1]:
        raise ValueError(f"Gold data cover {wavelengths[0]:.4f}–{wavelengths[-1]:.4f} µm only")
    return complex(np.interp(wavelength_um, wavelengths, n),
                   np.interp(wavelength_um, wavelengths, k))


def silica_n(wavelength_um: float) -> complex:
    """Malitson fused-silica Sellmeier equation; 20 °C, 0.21–6.7 µm."""
    if not 0.21 <= wavelength_um <= 6.7:
        raise ValueError("Fused-silica Malitson formula covers 0.21–6.7 µm")
    l2 = wavelength_um**2
    n2 = (1 + .6961663*l2/(l2-.0684043**2)
          + .4079426*l2/(l2-.1162414**2)
          + .8974794*l2/(l2-9.896161**2))
    return complex(np.sqrt(n2))


def bk7_n(wavelength_um: float) -> complex:
    """SCHOTT N-BK7 Sellmeier dispersion at 20 °C, 0.3–2.5 µm."""
    if not 0.3 <= wavelength_um <= 2.5:
        raise ValueError("N-BK7 formula covers 0.3–2.5 µm")
    l2 = wavelength_um**2
    n2 = (1 + 1.03961212*l2/(l2-.00600069867)
          + .231792344*l2/(l2-.0200179144)
          + 1.01046945*l2/(l2-103.560653))
    return complex(np.sqrt(n2))


def _schott_glass_n(wavelength_um: float, name: str,
                    terms: tuple[tuple[float, float], ...]) -> complex:
    """SCHOTT three-term Sellmeier form where C is already in µm²."""
    if not .365 <= wavelength_um <= 2.5:
        raise ValueError(f"{name} formula is limited here to 0.365–2.5 µm")
    l2 = wavelength_um*wavelength_um
    n2 = 1+sum(b*l2/(l2-c) for b, c in terms)
    if n2 <= 0 or not np.isfinite(n2):
        raise ValueError(f"{name} dispersion is invalid at this wavelength")
    return complex(np.sqrt(n2))


def caf2_n(wavelength_um: float) -> complex:
    """Malitson calcium-fluoride Sellmeier dispersion, 0.23–9.7 µm."""
    if not 0.23 <= wavelength_um <= 9.7:
        raise ValueError("CaF2 Malitson formula covers 0.23–9.7 µm")
    l2 = wavelength_um**2
    n2 = (1 + .5675888*l2/(l2-.050263605**2)
          + .4710914*l2/(l2-.1003909**2)
          + 3.8484723*l2/(l2-34.649040**2))
    return complex(np.sqrt(n2))


@lru_cache(maxsize=1)
def silicon_table() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    for line in SILICON_DATA.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            rows.append(tuple(float(part) for part in parts))
        except ValueError:
            continue
    array = np.asarray(rows)
    if array.shape[1] != 3 or len(array) < 20:
        raise RuntimeError("Silicon optical-constant table is incomplete")
    return array[:,0], array[:,1], array[:,2]


def silicon_n(wavelength_um: float) -> complex:
    """Green 2008 intrinsic-Si n,k at 300 K; no extrapolation."""
    wavelengths, n, k = silicon_table()
    if not wavelengths[0] <= wavelength_um <= wavelengths[-1]:
        raise ValueError(f"Green silicon data cover {wavelengths[0]:.2f}–{wavelengths[-1]:.2f} µm")
    return complex(np.interp(wavelength_um, wavelengths, n),
                   np.interp(wavelength_um, wavelengths, k))


def material_n(name: str, wavelength_um: float, constant_n: float = 1.0) -> complex:
    if name.startswith("custom:"):
        record = _user_materials().get(name)
        if record is None:
            raise ValueError("Unknown imported material")
        values = np.asarray(record["rows"], dtype=float)
        if not values[0, 0] <= wavelength_um <= values[-1, 0]:
            raise ValueError(f"{record['name']} data cover {values[0,0]}–{values[-1,0]} µm only")
        return complex(np.interp(wavelength_um, values[:,0], values[:,1]),
                       np.interp(wavelength_um, values[:,0], values[:,2]))
    if name == "dielectric":
        if not np.isfinite(constant_n) or constant_n <= 0:
            raise ValueError("Constant refractive index must be positive and finite")
        return complex(constant_n)
    if name == "air":
        return 1+0j
    if name == "gold":
        return gold_n(wavelength_um)
    if name == "silica":
        return silica_n(wavelength_um)
    if name == "bk7":
        return bk7_n(wavelength_um)
    if name == "n_sf11":
        return _schott_glass_n(wavelength_um, "N-SF11", (
            (1.73759695, .013188707), (.313747346, .0623068142),
            (1.89878101, 155.23629)))
    if name == "n_bak4":
        return _schott_glass_n(wavelength_um, "N-BAK4", (
            (1.28834642, .00779980626), (.132817724, .0315631177),
            (.945395373, 105.965875)))
    if name == "caf2":
        return caf2_n(wavelength_um)
    if name == "silicon":
        return silicon_n(wavelength_um)
    if name == "baf2":
        return _sellmeier(wavelength_um, .2652, 10.346, 0,
            ((.643356,.057789),(.506762,.10968),(3.8261,46.3864)), "BaF2")
    if name == "lif":
        return _sellmeier(wavelength_um, .1, 11, 0,
            ((.92549,.07376),(6.96747,32.790)), "LiF")
    if name == "kf":
        return _sellmeier(wavelength_um, .15, 22, .55083,
            ((.29162,.126),(3.60001,51.55)), "KF")
    if name == "nacl":
        return _sellmeier(wavelength_um, .2, 30, .00055,
            ((.198,.050),(.48398,.100),(.38696,.128),(.25998,.158),
             (.08796,40.50),(3.17064,60.98),(.30038,120.34)), "NaCl")
    if name == "tio2":
        if not .43 <= wavelength_um <= 1.53:
            raise ValueError("Rutile TiO2 DeVore formula covers 0.43–1.53 µm")
        square = wavelength_um * wavelength_um
        value = 5.913 + .2441 / (square - .0803)
        if value <= 0 or not np.isfinite(value):
            raise ValueError("Rutile TiO2 dispersion is invalid at this wavelength")
        return complex(np.sqrt(value))
    if name == "al2o3":
        return _sellmeier(wavelength_um, .2, 5.0, 0,
            ((1.4313493,.0726631),(.65054713,.1193242),(5.3414021,18.028251)),
            "alpha-Al2O3 sapphire ordinary ray")
    if name == "mgf2":
        return _sellmeier(wavelength_um, .2, 7.0, 0,
            ((.48755108,.04338408),(.39875031,.09461442),(2.3120353,23.793604)),
            "MgF2 ordinary ray")
    if name == "sin":
        if not .45 <= wavelength_um <= 1.65:
            raise ValueError("Silicon-nitride Sellmeier fit is limited here to 0.45–1.65 µm")
        square = wavelength_um*wavelength_um
        return complex(np.sqrt(1+2.9277*square/(square-.1321**2)
                               +.0200*square/(square-1.8703**2)))
    if name in RAKIC_LD:
        return _rakic_n(name, wavelength_um)
    if name in POLYMER_INDICES:
        if not POLYMER_WAVELENGTHS[0] <= wavelength_um <= POLYMER_WAVELENGTHS[-1]:
            raise ValueError("Sultanova polymer measurements cover 0.4368–1.052 µm")
        return complex(np.interp(wavelength_um, POLYMER_WAVELENGTHS,
                                 POLYMER_INDICES[name]))
    raise ValueError("Unknown material")


def material_energy_terms(name: str, wavelength_um: float,
                          constant_n: float = 1.0) -> dict:
    """Return local epsilon terms used for loss and weak-loss stored energy.

    For nonmagnetic isotropic media, d(omega*epsilon')/domega equals
    epsilon' - lambda*d(epsilon')/dlambda. The stored-energy expression is
    only labelled valid when loss is small and the derivative is positive.
    """
    center = material_n(name, wavelength_um, constant_n)**2
    step = max(1e-6, abs(wavelength_um)*1e-3)
    derivative = None
    for _ in range(12):
        try:
            lower = material_n(name, wavelength_um-step, constant_n)**2
            upper = material_n(name, wavelength_um+step, constant_n)**2
            derivative = (upper.real-lower.real)/(2*step)
            break
        except ValueError:
            step /= 2
    if derivative is None:
        raise ValueError("Material range is too narrow for an energy-density derivative")
    coefficient = float(center.real-wavelength_um*derivative)
    index = material_n(name, wavelength_um, constant_n)
    loss_ratio = float(abs(index.imag)/max(abs(index.real), 1e-15))
    valid = bool(coefficient > 0 and loss_ratio <= .01)
    return {"epsilon_real": float(center.real), "epsilon_imag": float(center.imag),
            "electric_energy_coefficient": coefficient,
            "derivative_step_um": step, "loss_ratio_k_over_n": loss_ratio,
            "weak_loss_valid": valid}


def material_catalog(wavelength_um: float) -> dict:
    """Optical properties and provenance for the novice-facing material inspector."""
    if not np.isfinite(wavelength_um) or wavelength_um <= 0:
        raise ValueError("Inspection wavelength must be finite and positive")
    items = []
    custom = {key: (item["name"], f"{item['rows'][0][0]}–{item['rows'][-1][0]} µm",
                    item["source"]) for key, item in _user_materials().items()}
    for key, (label, range_label, source) in {**MATERIAL_INFO, **custom}.items():
        if key.startswith("custom:"):
            category = "Imported"
        elif key in {"gold", "silver_rakic", "aluminum_rakic", "copper_rakic"}:
            category = "Metals"
        elif key in POLYMER_INDICES:
            category = "Polymers"
        elif key in {"silicon", "sin", "tio2", "al2o3"}:
            category = "Semiconductors and oxides"
        elif key in {"silica", "bk7", "n_sf11", "n_bak4", "caf2", "baf2", "lif", "kf", "nacl", "mgf2"}:
            category = "Glasses and crystals"
        else:
            category = "Reference media"
        row = {"key": key, "name": label, "category": category,
               "range": range_label, "source": source}
        try:
            value = material_n(key, wavelength_um)
            row.update({"n": float(value.real), "k": float(value.imag), "available": True})
        except ValueError:
            row.update({"n": None, "k": None, "available": False})
        items.append(row)
    return {"wavelength_um": wavelength_um, "materials": items,
            "note": "n+ik is the complex refractive index. k controls attenuation; Sellmeier entries here omit absorption."}
