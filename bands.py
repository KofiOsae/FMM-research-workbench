"""2D scalar plane-wave expansion for infinite, lossless, nonmagnetic media.

TM means E_z, TE means H_z (the MPB convention). This is a 2D infinite-crystal model; its bands
are not the resonances of a finite-thickness patterned slab.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import replace

from run_jobs import progress

import numpy as np
from scipy.linalg import eigh


@dataclass(frozen=True)
class BandModel:
    period_um: float = 0.6
    lattice: str = "square"  # square, rectangular, hexagonal (triangular Bravais)
    aspect_ratio: float = 1.0  # |a2| / |a1| for rectangular only
    background_n: float = 1.0
    inclusion_n: float = 3.0
    shape: str = "disk"
    fill: float = 0.25  # radius / period for disk; side / period for square
    ring_inner: float = 0.14
    fourier_order: int = 3
    grid_size: int = 256
    points_per_segment: int = 20
    bands: int = 8

    def validate(self) -> None:
        values = (self.period_um, self.aspect_ratio, self.background_n, self.inclusion_n, self.fill, self.ring_inner)
        if not all(np.isfinite(values)) or min(self.period_um, self.aspect_ratio,
                                                self.background_n, self.inclusion_n) <= 0:
            raise ValueError("Period and refractive indices must be finite and positive")
        if self.shape not in ("disk", "square", "ring"):
            raise ValueError("Shape must be disk, square, or ring")
        if self.lattice not in ("square", "rectangular", "hexagonal"):
            raise ValueError("Lattice must be square, rectangular, or hexagonal")
        if not .25 <= self.aspect_ratio <= 4:
            raise ValueError("Rectangular lattice aspect ratio must be 0.25–4")
        upper = 1 if self.shape == "square" else .5
        if not 0 < self.fill < upper or (self.shape == "ring" and not 0 < self.ring_inner < self.fill):
            raise ValueError("Square side must be below one period; disk/ring outer radius below half a period")
        if not 1 <= self.fourier_order <= 15:
            raise ValueError("Fourier order must be 1–15")
        if not 32 <= self.grid_size <= 1024:
            raise ValueError("Grid size must be 32–1024")
        if not 3 <= self.points_per_segment <= 100:
            raise ValueError("Points per segment must be 3–100")
        if not 1 <= self.bands <= (2*self.fourier_order+1)**2:
            raise ValueError("Requested bands exceed the Fourier basis")


def dielectric_grid(model: BandModel) -> np.ndarray:
    model.validate()
    u = (np.arange(model.grid_size) + 0.5) / model.grid_size - 0.5
    uu, vv = np.meshgrid(u, u, indexing="ij")
    if model.lattice == "hexagonal":
        x, y = uu + .5*vv, np.sqrt(3)/2*vv
    elif model.lattice == "rectangular":
        x, y = uu, model.aspect_ratio*vv
    else:
        x, y = uu, vv
    r = np.sqrt(x*x + y*y)
    if model.shape == "disk":
        inside = r < model.fill
    elif model.shape == "square":
        inside = (np.abs(x) < model.fill/2) & (np.abs(y) < model.fill/2)
    else:
        inside = (r < model.fill) & (r > model.ring_inner)
    return np.where(inside, model.inclusion_n**2, model.background_n**2)


def _convolution(values: np.ndarray, mx: np.ndarray, my: np.ndarray) -> np.ndarray:
    fft = np.fft.fft2(values) / values.size
    return fft[(mx[:, None]-mx[None, :]) % values.shape[0],
               (my[:, None]-my[None, :]) % values.shape[1]]


def _lattice(model: BandModel) -> tuple[np.ndarray, np.ndarray, list[str]]:
    if model.lattice == "hexagonal":
        direct = np.array([[1, .5], [0, np.sqrt(3)/2]])
        vertices, labels = np.array([[0,0],[.5,0],[2/3,1/3],[0,0.]]), ["Γ","M","K","Γ"]
    else:
        ratio = model.aspect_ratio if model.lattice == "rectangular" else 1
        direct = np.array([[1,0],[0,ratio]])
        vertices, labels = np.array([[0,0],[.5,0],[.5,.5],[0,0.]]), ["Γ","X","M","Γ"]
    return np.linalg.inv(direct).T, vertices, labels


def _path(model: BandModel) -> tuple[np.ndarray, np.ndarray, list[int], list[str]]:
    reciprocal, vertices, labels = _lattice(model)
    parts = [np.linspace(vertices[i], vertices[i+1], model.points_per_segment, endpoint=False)
             for i in range(3)]
    points = np.vstack((*parts, vertices[-1:]))
    cartesian = points @ reciprocal.T
    distances = np.r_[0, np.cumsum(np.linalg.norm(np.diff(cartesian, axis=0), axis=1))]
    return points, distances, [i*model.points_per_segment for i in range(4)], labels


def solve_bands(model: BandModel) -> dict:
    """Return normalized frequency a/lambda on Γ-X-M-Γ for TE and TM."""
    eps = dielectric_grid(model)
    orders = np.arange(-model.fourier_order, model.fourier_order+1)
    mx, my = np.meshgrid(orders, orders, indexing="ij")
    mx, my = mx.ravel(), my.ravel()
    eps_conv = _convolution(eps, mx, my)
    inverse_conv = _convolution(1/eps, mx, my)
    eps_conv = (eps_conv + eps_conv.conj().T)/2
    inverse_conv = (inverse_conv + inverse_conv.conj().T)/2
    reciprocal, _, _ = _lattice(model)
    points, distances, ticks, labels = _path(model)
    te, tm = [], []
    for point_index, (kx, ky) in enumerate(points):
        progress(point_index, len(points), "Bloch wavevectors")
        vectors = np.column_stack((mx+kx, my+ky)) @ reciprocal.T
        gx, gy = vectors[:,0], vectors[:,1]
        diagonal = np.diag(gx*gx+gy*gy)
        tm_values = eigh(diagonal, eps_conv, eigvals_only=True, subset_by_index=(0,model.bands-1))
        tm_matrix = gx[:,None]*inverse_conv*gx[None,:] + gy[:,None]*inverse_conv*gy[None,:]
        te_values = eigh(tm_matrix, eigvals_only=True, subset_by_index=(0,model.bands-1))
        te.append(np.sqrt(np.maximum(te_values, 0)).tolist())
        tm.append(np.sqrt(np.maximum(tm_values, 0)).tolist())
    return {
        "distance": distances.tolist(), "ticks": ticks,
        "tick_labels": labels,
        "TE": te, "TM": tm,
        "normalized_frequency": "a / vacuum wavelength",
        "unit_cell_epsilon": eps[::max(1,model.grid_size//64), ::max(1,model.grid_size//64)].tolist(),
    }


def full_zone_gaps(model: BandModel, points: int = 9) -> dict:
    """Screen complete TE/TM gaps over one reciprocal primitive cell."""
    model.validate()
    if not 5 <= points <= 31 or points % 2 == 0:
        raise ValueError("Use an odd full-zone grid from 5 to 31 points per reciprocal axis")
    if model.fourier_order > 7:
        raise ValueError("Interactive full-zone screening is capped at Fourier order 7")
    eps = dielectric_grid(model)
    orders = np.arange(-model.fourier_order, model.fourier_order+1)
    mx, my = np.meshgrid(orders, orders, indexing="ij"); mx, my = mx.ravel(), my.ravel()
    eps_conv = _convolution(eps, mx, my); inverse_conv = _convolution(1/eps, mx, my)
    eps_conv = (eps_conv+eps_conv.conj().T)/2
    inverse_conv = (inverse_conv+inverse_conv.conj().T)/2
    reciprocal, _, _ = _lattice(model)
    axis = np.linspace(-.5, .5, points)
    values = {"TE": [], "TM": []}
    total = points*points
    for index, (kx, ky) in enumerate((x, y) for y in axis for x in axis):
        progress(index, total, "Full Brillouin-zone screening")
        vectors = np.column_stack((mx+kx, my+ky)) @ reciprocal.T
        gx, gy = vectors[:,0], vectors[:,1]
        diagonal = np.diag(gx*gx+gy*gy)
        tm = eigh(diagonal, eps_conv, eigvals_only=True, subset_by_index=(0,model.bands-1))
        matrix = gx[:,None]*inverse_conv*gx[None,:]+gy[:,None]*inverse_conv*gy[None,:]
        te = eigh(matrix, eigvals_only=True, subset_by_index=(0,model.bands-1))
        values["TE"].append(np.sqrt(np.maximum(te, 0)))
        values["TM"].append(np.sqrt(np.maximum(tm, 0)))
    extrema, gaps = {}, []
    for polarization in ("TE", "TM"):
        array = np.asarray(values[polarization]).reshape(points, points, model.bands)
        minimum, maximum = array.min(axis=(0,1)), array.max(axis=(0,1))
        extrema[polarization] = {"minimum":minimum.tolist(), "maximum":maximum.tolist()}
        for band in range(model.bands-1):
            lower, upper = float(maximum[band]), float(minimum[band+1])
            if upper > lower:
                gaps.append({"polarization":polarization,"between_bands":[band+1,band+2],
                             "lower_a_over_lambda":lower,"upper_a_over_lambda":upper,
                             "relative_width":2*(upper-lower)/(upper+lower),
                             "wavelength_um":[model.period_um/upper, model.period_um/lower]})
    return {"grid_points_per_axis":points,"reduced_k_axis":axis.tolist(),
            "sampled_reciprocal_cell":"-0.5 ≤ k1,k2 ≤ 0.5; boundary duplicates are retained for audit",
            "band_extrema":extrema,"complete_sampled_gaps":gaps,
            "warning":"A sampled complete-gap screen. Repeat with a denser k grid, higher plane-wave order, and material uncertainty before making a complete-gap claim."}


def band_convergence(model: BandModel, base_result: dict | None = None,
                     tolerance: float = .01) -> dict:
    base = base_result if base_result is not None else solve_bands(model)
    refined = solve_bands(replace(model, fourier_order=model.fourier_order+1))
    changes = {key: np.max(np.abs(np.asarray(base[key])-np.asarray(refined[key])), axis=0).tolist()
               for key in ("TE", "TM")}
    return {"order_pair": [model.fourier_order, model.fourier_order+1],
            "max_change_by_band": changes, "tolerance": tolerance,
            "converged_by_band": {key: [delta <= tolerance for delta in values]
                                  for key, values in changes.items()}}


def solve_band_mode(model: BandModel, polarization: str, kx: float, ky: float,
                    band: int, resolution: int = 101) -> dict:
    """Reconstruct scalar Ez (TM) or Hz (TE) in one primitive cell."""
    model.validate()
    if polarization not in ("TE", "TM"):
        raise ValueError("Band mode polarization must be TE or TM")
    if not 0 <= band < model.bands:
        raise ValueError("Selected band is outside the requested band count")
    if not np.isfinite([kx, ky]).all() or not 33 <= resolution <= 257:
        raise ValueError("Use finite reduced k coordinates and 33–257 field samples")
    eps = dielectric_grid(model)
    orders = np.arange(-model.fourier_order, model.fourier_order+1)
    mx, my = np.meshgrid(orders, orders, indexing="ij")
    mx, my = mx.ravel(), my.ravel()
    eps_conv = _convolution(eps, mx, my)
    inverse_conv = _convolution(1/eps, mx, my)
    eps_conv = (eps_conv+eps_conv.conj().T)/2
    inverse_conv = (inverse_conv+inverse_conv.conj().T)/2
    reciprocal, _, _ = _lattice(model)
    vectors = np.column_stack((mx+kx, my+ky)) @ reciprocal.T
    gx, gy = vectors[:, 0], vectors[:, 1]
    if polarization == "TM":
        matrix = np.diag(gx*gx+gy*gy)
        values, coefficients = eigh(matrix, eps_conv,
                                    subset_by_index=(0, model.bands-1))
    else:
        matrix = gx[:, None]*inverse_conv*gx[None, :] + gy[:, None]*inverse_conv*gy[None, :]
        values, coefficients = eigh(matrix, subset_by_index=(0, model.bands-1))
    coefficient = coefficients[:, band]
    axis = np.linspace(-.5, .5, resolution, endpoint=False)
    uu, vv = np.meshgrid(axis, axis, indexing="ij")
    phase = np.exp(2j*np.pi*((mx+kx)[:, None, None]*uu[None, :, :]
                             +(my+ky)[:, None, None]*vv[None, :, :]))
    field = np.tensordot(coefficient, phase, axes=(0, 0))
    field /= max(float(np.max(np.abs(field))), 1e-30)
    intensity = np.abs(field)**2
    inversion_overlap = None
    trim = np.allclose(2*np.asarray([kx, ky]), np.round(2*np.asarray([kx, ky])), atol=1e-10)
    if trim:
        inverted = np.roll(np.roll(field[::-1, ::-1], 1, axis=0), 1, axis=1)
        inversion_overlap = complex(np.vdot(field, inverted)/np.vdot(field, field))
    return {"polarization": polarization, "component": "H_z" if polarization == "TE" else "E_z",
            "band": band, "k_reduced": [float(kx), float(ky)],
            "normalized_frequency": float(np.sqrt(max(values[band], 0))),
            "u": axis.tolist(), "v": axis.tolist(),
            "real": field.real.tolist(), "imag": field.imag.tolist(),
            "intensity": intensity.tolist(), "phase": np.angle(field).tolist(),
            "epsilon": dielectric_grid(replace(model, grid_size=resolution)).tolist(),
            "inversion_overlap": None if inversion_overlap is None else
                {"real": float(inversion_overlap.real), "imag": float(inversion_overlap.imag),
                 "magnitude": float(abs(inversion_overlap))},
            "note": "The scalar field is periodic Bloch amplitude times its reduced-k phase in one primitive cell. Overall complex phase is arbitrary."}
