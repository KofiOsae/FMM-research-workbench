"""Angle–wavelength and incident-k-space maps for the finite FMM stack."""

from dataclasses import asdict
import math

import numpy as np

from stack import StackModel, solve_stack, stack_convergence


QUANTITIES = ("R", "T", "A", "R0", "T0")
POLARIZATION_QUANTITIES = ("power", "Ep_real", "Ep_imag", "Es_real", "Es_imag",
                           "S1", "S2", "S3", "orientation_deg", "ellipticity_deg")


def _validated_quantity(quantity: str) -> str:
    if quantity not in QUANTITIES:
        raise ValueError("Map quantity must be R, T, A, R0, or T0")
    return quantity


def angle_wavelength_map(model: StackModel, wavelength_start: float,
                         wavelength_stop: float, wavelength_points: int,
                         theta_start: float, theta_stop: float,
                         theta_points: int, quantity: str = "R") -> dict:
    """Scan vacuum wavelength and polar incidence angle at fixed azimuth."""
    quantity = _validated_quantity(quantity)
    values = np.asarray([wavelength_start, wavelength_stop, theta_start, theta_stop], float)
    if not np.isfinite(values).all() or wavelength_start <= 0 or wavelength_start >= wavelength_stop:
        raise ValueError("Use positive, increasing finite wavelength limits")
    if not -89 < theta_start < theta_stop < 89:
        raise ValueError("Use increasing signed incidence angles between -89° and +89°")
    if not 3 <= wavelength_points <= 41 or not 3 <= theta_points <= 41:
        raise ValueError("Use 3–41 points along each map axis")
    wavelengths = np.linspace(wavelength_start, wavelength_stop, wavelength_points)
    angles = np.linspace(theta_start, theta_stop, theta_points)
    grid = np.full((theta_points, wavelength_points), np.nan)
    failures = []
    for iy, theta in enumerate(angles):
        for ix, wavelength in enumerate(wavelengths):
            # A negative signed angle is the equivalent positive polar angle
            # with the in-plane propagation direction reversed by 180 degrees.
            signed_phi = model.phi_deg if theta >= 0 else model.phi_deg + 180
            current = StackModel(**{**asdict(model), "layers": model.layers,
                                    "wavelength_um": float(wavelength),
                                    "theta_deg": float(abs(theta)), "phi_deg": float(signed_phi)})
            try:
                grid[iy, ix] = solve_stack(current)[quantity]
            except (ValueError, np.linalg.LinAlgError) as exc:
                failures.append({"theta_deg": float(theta), "wavelength_um": float(wavelength),
                                 "reason": str(exc)})
    center = StackModel(**{**asdict(model), "layers": model.layers,
                           "wavelength_um": float(wavelengths[len(wavelengths)//2]),
                           "theta_deg": float(abs(angles[len(angles)//2])),
                           "phi_deg": float(model.phi_deg if angles[len(angles)//2] >= 0 else model.phi_deg+180)})
    convergence = stack_convergence(center)
    return {"kind": "angle_wavelength", "quantity": quantity,
            "wavelength_um": wavelengths.tolist(), "theta_deg": angles.tolist(),
            "values": [[None if not np.isfinite(v) else float(v) for v in row] for row in grid],
            "failed_points": failures, "center_convergence": convergence,
            "coordinates": "signed theta is measured across the surface normal; negative theta reverses the in-plane wavevector at phi+180 degrees"}


def kspace_map(model: StackModel, wavelength_um: float, rho_max: float,
               points: int, quantity: str = "R") -> dict:
    """Scan incident transverse wavevector at one vacuum wavelength.

    Coordinates u=kx/(n_inc*k0), v=ky/(n_inc*k0), so rho=sin(theta).
    Points outside the requested circular numerical aperture are null.
    """
    quantity = _validated_quantity(quantity)
    if not np.isfinite([wavelength_um, rho_max]).all() or wavelength_um <= 0:
        raise ValueError("Wavelength and k-space radius must be finite and positive")
    if not 0 < rho_max < math.sin(math.radians(89)):
        raise ValueError("k-space radius must satisfy 0 < k_parallel/(n_inc k0) < sin(89°)")
    if not 5 <= points <= 41 or points % 2 == 0:
        raise ValueError("Use an odd k-space grid from 5 to 41 points")
    axis = np.linspace(-rho_max, rho_max, points)
    grid = np.full((points, points), np.nan)
    failures = []
    for iy, v in enumerate(axis):
        for ix, u in enumerate(axis):
            rho = math.hypot(u, v)
            if rho > rho_max:
                continue
            theta = math.degrees(math.asin(rho))
            phi = 0.0 if rho < 1e-14 else math.degrees(math.atan2(v, u))
            current = StackModel(**{**asdict(model), "layers": model.layers,
                                    "wavelength_um": float(wavelength_um),
                                    "theta_deg": theta, "phi_deg": phi})
            try:
                grid[iy, ix] = solve_stack(current)[quantity]
            except (ValueError, np.linalg.LinAlgError) as exc:
                failures.append({"u": float(u), "v": float(v), "reason": str(exc)})
    center = StackModel(**{**asdict(model), "layers": model.layers,
                           "wavelength_um": float(wavelength_um),
                           "theta_deg": 0.0, "phi_deg": 0.0})
    convergence = stack_convergence(center)
    return {"kind": "kspace", "quantity": quantity, "u": axis.tolist(), "v": axis.tolist(),
            "values": [[None if not np.isfinite(value) else float(value) for value in row]
                       for row in grid], "failed_points": failures,
            "center_convergence": convergence,
            "coordinates": "u=kx/(n_inc k0), v=ky/(n_inc k0), rho=sin(theta); s/p rotate with phi"}


def polarization_kspace_map(model: StackModel, wavelength_um: float, rho_max: float,
                            points: int, port: str = "reflected",
                            order_m: int = 0, order_n: int = 0,
                            min_power: float = 1e-10) -> dict:
    """Map one coherent outgoing diffraction channel over incident k-space.

    Jones coefficients and Stokes parameters use each outgoing order's local
    right-handed p/s basis.  This is a driven-scattering observable, not the
    radiation coefficient of an isolated slab eigenmode.
    """
    if model.polarization not in ("s", "p"):
        raise ValueError("Far-field Jones/Stokes maps require coherent s or p input")
    if port not in ("reflected", "transmitted"):
        raise ValueError("Far-field port must be reflected or transmitted")
    if not np.isfinite([wavelength_um, rho_max, min_power]).all() or wavelength_um <= 0:
        raise ValueError("Wavelength, k-space radius, and power threshold must be finite")
    if not 0 < rho_max < math.sin(math.radians(89)):
        raise ValueError("k-space radius must satisfy 0 < k_parallel/(n_inc k0) < sin(89°)")
    if not 5 <= points <= 21 or points % 2 == 0:
        raise ValueError("Use an odd polarization grid from 5 to 21 points")
    if not 0 <= min_power < 1:
        raise ValueError("Minimum channel power must satisfy 0 <= threshold < 1")
    if abs(order_m) > 20 or abs(order_n) > 20:
        raise ValueError("Diffraction order indices must be between -20 and 20")

    axis = np.linspace(-rho_max, rho_max, points)
    maps = {key: np.full((points, points), np.nan) for key in POLARIZATION_QUANTITIES}
    valid = np.zeros((points, points), dtype=bool)
    failures = []
    power_key = "R" if port == "reflected" else "T"
    pol_key = port + "_polarization"
    for iy, v in enumerate(axis):
        for ix, u in enumerate(axis):
            rho = math.hypot(u, v)
            if rho > rho_max:
                continue
            theta = math.degrees(math.asin(rho))
            phi = 0.0 if rho < 1e-14 else math.degrees(math.atan2(v, u))
            current = StackModel(**{**asdict(model), "layers": model.layers,
                                    "wavelength_um": float(wavelength_um),
                                    "theta_deg": theta, "phi_deg": phi})
            try:
                result = solve_stack(current)
                channel = next((row for row in result["order_amplitudes"]
                                if row["m"] == order_m and row["n"] == order_n), None)
                if channel is None or channel[power_key] < min_power:
                    continue
                pol = channel[pol_key]
                maps["power"][iy, ix] = channel[power_key]
                scale = math.sqrt(max(pol["S0"], 1e-30))
                for component in ("Ep", "Es"):
                    maps[component + "_real"][iy, ix] = pol[component]["real"]/scale
                    maps[component + "_imag"][iy, ix] = pol[component]["imag"]/scale
                for key in ("S1", "S2", "S3"):
                    maps[key][iy, ix] = pol["normalized"][key]
                maps["orientation_deg"][iy, ix] = pol["orientation_deg"]
                maps["ellipticity_deg"][iy, ix] = pol["ellipticity_deg"]
                valid[iy, ix] = True
            except (ValueError, np.linalg.LinAlgError) as exc:
                failures.append({"u": float(u), "v": float(v), "reason": str(exc)})

    center = StackModel(**{**asdict(model), "layers": model.layers,
                           "wavelength_um": float(wavelength_um),
                           "theta_deg": 0.0, "phi_deg": 0.0})
    convergence = stack_convergence(center)
    serial = lambda array: [[None if not np.isfinite(value) else float(value)
                             for value in row] for row in array]
    return {"kind": "polarization_kspace", "u": axis.tolist(), "v": axis.tolist(),
            "maps": {key: serial(value) for key, value in maps.items()},
            "valid": valid.tolist(), "valid_points": int(valid.sum()),
            "total_inside_aperture": int(sum(math.hypot(u, v) <= rho_max
                                              for u in axis for v in axis)),
            "port": port, "order": {"m": int(order_m), "n": int(order_n)},
            "wavelength_um": float(wavelength_um), "min_power": float(min_power),
            "basis": "local right-handed p,s for each outgoing diffraction order",
            "jones_normalization": "Ep and Es are divided by sqrt(S0); phase retains the incident-wave reference",
            "stokes_convention": "S3=-2 Im(Ep Es*)",
            "failed_points": failures, "center_convergence": convergence,
            "warning": "Driven-scattering polarization includes the incident/background response; it is not by itself an eigenmode radiation vortex or a BIC topological charge."}
